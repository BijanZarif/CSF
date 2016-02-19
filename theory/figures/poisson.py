from dolfin import *

mesh = UnitIntervalMesh(10)

V = FunctionSpace(mesh,'CG',2)

u = TrialFunction(V)
v = TestFunction(V)

u0 = Constant(0)
u1 = Constant(1)
g = Constant(3.0)
f = Expression('x[0]')


class boundary0(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[0],0.0)
class boundary1(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[0],1.0)

SD = FacetFunction('size_t',mesh)
boundary0().mark(SD,0)
boundary1().mark(SD,1)

bc0 = DirichletBC(V,u0,SD,0)
bc1 = DirichletBC(V,u1,SD,1)
bcs = [bc0,bc1]


F = -inner(grad(u),grad(v))*dx + inner(f,v)*dx# + inner(Constant(2),v)*ds(1)

u = Function(V)

solve(lhs(F) == rhs(F), u, bcs)
plot(u,interactive=True)


