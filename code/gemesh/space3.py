from dolfin import *
N = 20
mesh = UnitSquareMesh(N,N)

class Boundary(SubDomain):
    def inside(self,x,on_bnd):
        return on_bnd

class Boundary_L(SubDomain):
    def inside(self,x,on_bnd):
        return near(x[0],0.0) and on_bnd

class Boundary_R(SubDomain):
    def inside(self,x,on_bnd):
        return near(x[0],1.0) and on_bnd

class Walls(SubDomain):
    def inside(self,x,on_bnd):
        return (near(x[1],0.0) or near(x[1],1.0)) and on_bnd

b = Boundary()
bl = Boundary_L()
br = Boundary_R()
wall = Walls()



bnd = FacetFunction('size_t',mesh)
bnd.set_all(0)
bl.mark(bnd,1)
br.mark(bnd,2)
wall.mark(bnd,3)

ds = Measure('ds')[bnd]

N = 20

V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)

VQR = MixedFunctionSpace([V,Q,V])

upw = TrialFunction(VQR)
v,q,r = TestFunctions(VQR)


u,p,w = split(upw)

p1 = Constant(1.0)
u0 = Expression(('0.25*x[1]*(1-x[1])','0'))
w0 = Expression(('sin(x[0])','0.0'))
f = Constant((1.0,0.0))
n = FacetNormal(mesh)

# Conservation of:

# Momentum 
F_m = inner(grad(u),grad(v))*dx + inner(div(v),p)*dx - inner(v,p1*n)*ds(2)

# Mass
F_c = inner(div(u),q)*dx

# Mesh:
F_w = inner(grad(w),grad(r))*dx - inner(f,r)*dx

noslip = Constant((0.0,0.0))
bcu = DirichletBC(VQR.sub(0),u0,bl)
bcw = DirichletBC(VQR.sub(2),noslip,b)

bcs = [bcu,bcw]

UPW_ = Function(VQR)

F = F_m + F_c + F_w

solve(lhs(F)==rhs(F), UPW_, bcs)

u_,p_,w_ = split(UPW_)
print w_
plot(u_)
plot(p_)
plot(w_)
interactive()
