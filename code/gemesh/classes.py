from domains_and_boundaries import *
import numpy as np
parameters['allow_extrapolation'] = True


bot_out = bottom_outer()
bot_in = bottom_inner()
out_wall = outer_wall()
in_wall = inner_wall()
top_out = top_outer()
top_in = top_inner()
csc_bnd = csc_boundary()
in_csc = inside_csc()
cord = cord()
CSF_s = CSF_space()

CSF_mesh = SubMesh(mesh,CSF_s)

boundaries = FacetFunction('size_t',CSF_mesh) # CSF_mesh
boundaries.set_all(0)


out_wall.mark(boundaries,1)
in_wall.mark(boundaries,2)
top_out.mark(boundaries,3)
bot_out.mark(boundaries,4)
bot_in.mark(boundaries,5)
top_in.mark(boundaries,6)
csc_bnd.mark(boundaries,7)

'''
plot(boundaries)
interactive()
import sys
sys.exit()
'''
noslip = Constant((0.0,0.0))

pb = Constant(0.0)
P = 1000
pt = Expression(('P*(sin(2*pi*t) + 0.2*sin(8*pi*t))'),t=0,P=P)

flow = Expression(('0.0','(x[0]-x0)*(x[0]+x0)*(x[0]-x1)*(x[0]+x1)*pow(10,10)'),x0=x0,x1=x1)



V = VectorFunctionSpace(CSF_mesh,'CG',2)
Q = FunctionSpace(CSF_mesh,'CG',1)
W = TensorFunctionSpace(CSF_mesh,'CG',2)

u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)


inflow = DirichletBC(V,flow,boundaries, 3)
bcu0 = DirichletBC(V,noslip,boundaries, 1)
bcu1 = DirichletBC(V,noslip,boundaries, 2)

bcU0 = DirichletBC(V,noslip,boundaries,4)
bcU1 = DirichletBC(V,noslip,boundaries,3)

bcp1 = DirichletBC(Q,pt,boundaries, 3)
bcp2 = DirichletBC(Q, pb, boundaries ,4)


bcu = [bcu0,bcu1]#,bcu2,bcu3,bcu4]
bcp = [bcp1,bcp2]


ufile = File('results/vel.pvd')
pfile = File('results/pr.pvd')

t = 0
T = 0.4
dt = 0.001

u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx + inner(u,v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

disp = Function(V)

class MyEx(Expression):
    def __init__(self,mesh,tau):
        self.mesh = mesh
        self.tau = tau
    def eval_cell(self,values,x,c):
        tau = self.tau(x)
        cell = Cell(self.mesh,c.index)
        n = cell.normal(c.local_facet)
        values[0] = (tau[0]*n[0] + tau[1]*n[1])*0.001/P
        values[1] = (tau[2]*n[0] + tau[3]*n[1])*0.001/P
    def value_shape(self):
        return (2,)

while t < T:
    pt.t = t

    # Update pressure boundary condition

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", prec)
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")
    end()
	
    ufile << u1
    pfile << p1
    

    pI = Identity(2)*p1
    viscous = mu*(grad(u1)+(grad(u1)).T)

    tau = -pI + viscous
    tau = project(tau,W)
    d = MyEx(CSF_mesh,tau)
    bcwall = DirichletBC(V,d,boundaries,2)
    bcDisp = [bcu0,bcwall]
    F = inner(grad(u),grad(v))*dx + inner(Constant(('0.0','0.0')),v)*dx
    solve(lhs(F) == rhs(F), disp, bcs=bcDisp)

    CSF_mesh.move(disp)
    CSF_mesh.bounding_box_tree().build(CSF_mesh)
    plot(boundaries)
    t += dt
    print "t =", t
    u0.assign(u1)
File('initial_data/u_double_refine.xml') << u0

interactive()
