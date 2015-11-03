from dolfin import *

mesh = UnitIntervalMesh(5)

V = FunctionSpace(mesh,'CG',2)

u = TrialFunction(V)
v = TestFunction(V)

u0 = Constant(0)
u1 = Constant(1)
g = Constant(3.0)
f = Expression('2.0')
n = FacetNormal(mesh)

def boundary0(x,on_bnd):
	return on_bnd and near(x[0],0.0)
def boundary1(x,on_bnd):
	return on_bnd and near(x[0],1.0)

bc0 = DirichletBC(V,u0,boundary0)
bc1 = DirichletBC(V,u1,boundary1)
bcs = [bc0,bc1]

F = -inner(grad(u),grad(v))*dx + f*v*dx

u = Function(V)

solve(lhs(F) == rhs(F), u, bcs)
plot(u,interactive=True)


