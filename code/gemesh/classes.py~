from domains_and_boundaries import *


subdomains = MeshFunction('uint',mesh,2)

CSF_mesh = SubMesh(mesh,CSF_s)
Cord_mesh = SubMesh(mesh,cord)
CSC_mesh = SubMesh(mesh,in_csc)

CSF_s.mark(subdomains,0)
cord.mark(subdomains,1)
in_csc.mark(subdomains,0)

#boundaries = FacetFunction('size_t',mesh) # CSF_mesh
#boundaries.set_all(0)
boundaries = MeshFunction('uint',mesh ,mesh.topology().dim()-1)


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



V = VectorFunctionSpace(mesh,'CG',2)
U = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)
W = TensorFunctionSpace(mesh,'CG',2)

u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

d = TrialFunction(V)    # displacement function


inflow = DirichletBC(V,flow,boundaries, 3)  # top CSF
bcu0 = DirichletBC(V,noslip,boundaries, 1)  # outer wall
bcu1 = DirichletBC(V,noslip,boundaries, 2)  # inner wall
bcu2 = DirichletBC(V,noslip,boundaries, 5)
bcu3 = DirichletBC(V,noslip,boundaries, 6)
bcu4 = DirichletBC(V,noslip,boundaries, 7)

bcU0 = DirichletBC(V,noslip,boundaries,4)   # bot CSF
bcU1 = DirichletBC(V,noslip,boundaries,3)   # top CSF

bcp1 = DirichletBC(Q,pt,boundaries, 3)
bcp2 = DirichletBC(Q, pb, boundaries ,4)
bcp3 = DirichletBC(Q,pt, boundaries, 6)     # top cord
bcp4 = DirichletBC(Q,pb, boundaries, 5)     # bottom cord


bcd0 = DirichletBC(V,noslip,boundaries,1)
bcd1 = DirichletBC(V,noslip,boundaries,2)
bcd2 = DirichletBC(V,noslip,boundaries,3)
bcd3 = DirichletBC(V,noslip,boundaries,4)
bcd4 = DirichletBC(V,noslip,boundaries,5)
bcd5 = DirichletBC(V,noslip,boundaries,6)
bcd6 = DirichletBC(V,noslip,boundaries,1)

bcu = [bcu0,bcu1,bcu2,bcu3,bcu4]
bcp = [bcp1,bcp2,bcp3,bcp4]
bcd = [bcd0,bcd1,bcd2,bcd3,bcd4,bcd5,bcd6]


ufile = File('results/vel.pvd')
pfile = File('results/pr.pvd')
dfile = File('results/disp.pvd')

t = 0
T = 0.2
dt = 0.001

u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)
disp = Function(V)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

n = FacetNormal(mesh)

# Tentative velocity step 
F1 = (1/k)*inner(u - u0, v)*dx(0, subdomain_data=subdomains) + inner(grad(u0)*u0, v)*dx(0, subdomain_data=subdomains) + nu*inner(grad(u), grad(v))*dx(0, subdomain_data=subdomains) - inner(f, v)*dx(0, subdomain_data=subdomains) + inner(grad(d),grad(v))*dx(0,subdomain_data=subdomains)

lamda = 0.001
kappa = 1.0

#(mu+lamda)*inner(div(d),div(v))*dx(1,subdomain_data=subdomains) - (1.0/kappa)*inner(u,v)*dx(1,subdomain_data=subdomains) + inner(div(d),q)*dx(1,subdomain_data=subdomains) + dot(u,grad(q))*dx(1,subdomain_data=subdomains) - mu*inner(grad(u).T*n,v)*ds(3)

a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx(0, subdomain_data=subdomains) + inner(grad(p),grad(q))*dx(1, subdomain_data=subdomains)
L2 = -(1/k)*div(u1)*q*dx(0, subdomain_data=subdomains) + inner(Constant('0.1'),q)*dx(1, subdomain_data=subdomains)

# Velocity update
a3 = inner(u, v)*dx(0, subdomain_data=subdomains) + inner(u,v)*dx(0, subdomain_data=subdomains) + inner(grad(u),grad(v))*dx(1, subdomain_data=subdomains)
L3 = inner(u1, v)*dx(0, subdomain_data=subdomains) - k*inner(grad(p1), v)*dx(0, subdomain_data=subdomains) + inner(Constant(('0.1','0')),v)*dx(1, subdomain_data=subdomains)

a4 = mu*inner(grad(d),grad(v))*dx(1,subdomain_data=subdomains) + + mu*inner(grad(d),grad(v))*dx(0,subdomain_data=subdomains)

L4 = inner(Constant(('0.0','0.0')),v)*dx(1,subdomain_data=subdomains)+ inner(Constant('0.0','0.0')),v)*dx(0,subdomain_data=subdomains)
# Assemble matrices

A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)
A4 = assemble(a4)

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
        values[0] = (tau[0]*n[0] + tau[1]*n[1])*0.005/P
        values[1] = (tau[2]*n[0] + tau[3]*n[1])*0.005/P
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
	
	# Compute displacement
	begin("Computing displacement")
	b4 = assemble(L4)
	[bc.apply(A4,b4) for bc in bcd]
	solve(A4, disp.vector(),b4)
	end()
	
	dfile << disp
    ufile << u1
    pfile << p1
    
    '''
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
    '''
    t += dt
    print "t =", t
    u0.assign(u1)
#File('initial_data/u_double_refine.xml') << u0
#File('end_mesh.xml') << CSF_mesh

interactive()
