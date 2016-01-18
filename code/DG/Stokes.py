from dolfin import *
class Right(SubDomain):
	def inside(self, x, on_boundary):
		return abs(x[0] - 1.0) < DOLFIN_EPS and on_boundary

class Left(SubDomain):
	def inside(self, x, on_boundary):
		return abs(x[0]) < DOLFIN_EPS and on_boundary

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return abs(x[1] - 1.0) < DOLFIN_EPS and on_boundary

class Bottom(SubDomain):
	def inside(self, x, on_boundary):
		return abs(x[1]) < DOLFIN_EPS and on_boundary


# mesh
right = Right()
left = Left()
top = Top()
bottom = Bottom()

mesh = UnitSquareMesh(40,40)
boundaries = FacetFunction('size_t',mesh)
boundaries.set_all(0)
left.mark(boundaries,1)
right.mark(boundaries,2)

ds = Measure('ds')[boundaries]

# Defining the function spaces
Q = FunctionSpace(mesh, "CG", 1)
V = VectorFunctionSpace(mesh, "CG", 2)
VQ = V*Q
# Create velocity Function


# Test and trial functions
v,q  = TestFunctions(VQ)
up = TrialFunction(VQ)
u,p = split(up)



# Diffusivity
mu = Constant(0.1)

# Source term
f = Constant((0.0,0.0))

# Penalty term
alpha = Constant(5.0)

# Mesh-related functions
n = FacetNormal(mesh)
h = CellSize(mesh)
h_avg = (h('+') + h('-'))/2

# ( dot(v, n) + |dot(v, n)| )/2.0
un = (dot(u, n) + abs(dot(u, n)))/2.0

p_in = 10.0
p_out = 1.0

n = FacetNormal(mesh)

# Bilinear form
a_int = inner(grad(u),grad(v))*dx - inner(p,div(v))*dx - inner(q,div(u))*dx

'''
a_diff = mu('+')*(alpha('+')/h('+'))*dot(jump(v, n), jump(phi, n))*dS \
- mu('+')*dot(avg(grad(v)), jump(phi, n))*dS \
- mu('+')*dot(jump(v, n), avg(grad(phi)))*dS

a_p = dot(jump(v), un('+')*phi('+') - un('-')*phi('-') )*dS \
+ dot(v, un*phi)*ds
'''

a = a_int # + a_fac + a_vel


# Linear form
L = -inner(v,p_in*n)*ds(1) - inner(v,p_out*n)*ds(2)

A = a+L

# Set up boundary condition (apply strong BCs)
noslip = Constant((0.0,0.0))
u_analytical = Expression(('x[1]*(1-x[1])','0.0'))
p_analytical = Expression('-2+x[0]')

bcuT = DirichletBC(VQ.sub(0),noslip,top)
bcuB = DirichletBC(VQ.sub(0),noslip,bottom)
p_in = 10.0
p_out = 1.0


bcs = [bcuT, bcuB]
# Solution function
UP_ = Function(VQ)

# Assemble and apply boundary conditions
"""
A = assemble(a)
b = assemble(L)
[bc.apply(A, b) for bc in bcs]
"""
F = A+L

# Solve system

solve(lhs(F) ==	 rhs(F), UP_, bcs=bcs)
u,p = UP_.split(True)
plot(u)
plot(p)
interactive()

