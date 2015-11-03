from dolfin import *

mesh = UnitSquareMesh(10,10)

V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)

VQ = V*Q

up = TrialFunction(VQ)
v,q = TestFunctions(VQ)

noslip = Constant((0.0,0.0))

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0)
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)

left = Left()
right = Right()
bottom = Bottom()
top = Top()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)

left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
'''
plot(boundaries)
interactive()
import sys
sys.exit()
'''
bcu = DirichletBC(VQ.sub(0),noslip,left)
bcu2 = DirichletBC(VQ.sub(0),noslip,right)
bcp = DirichletBC(VQ.sub(1),Constant(100),left)
bcp2 = DirichletBC(VQ.sub(1),Constant(0),right)

pk = Expression('2000*x[0]')
ps = Expression('1000*x[1]')

mu = Constant(1E3)
lamda = Constant(1E5)
K = Constant(1E-5)

ds = Measure("ds")[boundaries]



n = FacetNormal(mesh)
u,p = split(up)

F_m = mu*inner(grad(u),grad(v))*dx + \
      (mu+lamda)*inner(div(u),div(v))*dx - \
      inner(p,div(v))*dx + \
      inner(v,pk*n)*ds(1)

F_c = K*inner(grad(p),grad(q))*dx - inner(q,K*ps)*ds(2)

F = F_m + F_c

UP = Function(VQ)

solve(lhs(F) == rhs(F), UP, bcs=[bcu,bcp])
u_,p_ = UP.split(True)
plot(u_)
plot(p_)
interactive()

