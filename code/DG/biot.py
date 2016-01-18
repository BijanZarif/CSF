from dolfin import *
import sys
set_log_active(False)

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

h0 = 1
e0 = 1
ep0 = 1

for N in [4,8,16,32]:#,64]:
	# mesh
	right = Right()
	left = Left()
	top = Top()
	bottom = Bottom()

	mesh = UnitSquareMesh(N,N)
	boundaries = FacetFunction('size_t',mesh)
	boundaries.set_all(0)
	left.mark(boundaries,1)
	right.mark(boundaries,2)
	top.mark(boundaries,3)
	bottom.mark(boundaries,3)

	ds = Measure('ds')[boundaries]

	# Defining the function spaces
	Q = FunctionSpace(mesh, "CG", 1)
	V = VectorFunctionSpace(mesh, "CG", 2)
	VQ = V*Q
	W = VectorFunctionSpace(mesh,'CG',3)
	X = FunctionSpace(mesh,'CG',3)
	# Create velocity Function


	# Test and trial functions
	v,q  = TestFunctions(VQ)
	u,p = TrialFunctions(VQ)
	# Diffusivity
	mu = Constant(1.0)
        lamda = Constant(1.1)
	# Source term
	f = Constant((0.0,0.0))

	# Penalty term
	alpha = Constant(14.0)
	beta = Constant(10.0)
        gamma = Constant(8.0)

	# Mesh-related functions
	n = FacetNormal(mesh)
	h = CellSize(mesh)
	h_avg = (h('+') + h('-'))/2

	# ( dot(v, n) + |dot(v, n)| )/2.0
	p_in = Expression('std::sin(x[1])')
	p_out = Expression('std::cos(x[1])')
        k = 1.0
	n = FacetNormal(mesh)

	def eps(u):
		return 0.5*(grad(u) + grad(u).T)

	# Bilinear form
     	a_int = mu*inner(grad(u),grad(v))*dx + inner(p,div(v))*dx + k*inner(grad(q),grad(p))*dx + (mu+lamda)*inner(div(u),div(v))*dx

	a_diff =  - mu('+')*inner(dot(avg(grad(u)),n('+')),jump(v))*dS \
                  - mu('+')*inner(dot(avg(grad(v)),n('+')),jump(u))*dS \
		  + mu('+')*beta('+')/h('+')*inner(jump(u),jump(v))*dS \
		  - inner(avg(p),dot(jump(v),n('+')))*dS \
                  - inner(avg(q),dot(jump(u),n('+')))*dS \
		  #+ mu*alpha/h*inner(u,v)*ds(3) #\
                  #- mu*inner(dot(grad(u),n),v)*ds(3) \
		  #- mu*inner(dot(grad(v),n),u)*ds(3) \
                  #- mu*inner(dot(grad(u),n),v)*ds(3) \
		  #- mu*inner(dot(grad(v),n),u)*ds(3) \

        a_biot =  - (mu('+')+lamda('+'))*inner(jump(v), avg(div(u))*n('+'))*dS \
                  - (mu('+')+lamda('+'))*inner(jump(u), avg(div(v))*n('+'))*dS \
                  + gamma('+')*(mu('+')+lamda('+'))/h('+')*inner(jump(u),jump(v))*dS \
                  - k*inner(avg(q),dot(jump(grad(p)),n('+')))*dS \
                  - k*inner(jump(q),dot(avg(grad(p)),n('+')))*dS 
        '''
                  - (mu+lamda)*inner(v, div(u)*n)*ds(1) \
                  - (mu+lamda)*inner(u, div(v)*n)*ds(1) \
                  - (mu+lamda)*inner(v, div(u)*n)*ds(2) \
                  - (mu+lamda)*inner(u, div(v)*n)*ds(2) \
                  + gamma*(mu+lamda)/h*inner(u,v)*ds(3)
        '''
	#a = a_int + a_diff# + a_fac + a_vel
        a = a_int + a_diff + a_biot

	# Linear form
        u_analytical = Expression(('std::sin(pi*x[0]) + std::cos(pi*x[1])','std::cos(pi*x[1]*x[0])'))
	p_analytical = Expression('std::sin(pi*x[0]*x[1])')
        U = project(u_analytical,W)
        P = project(p_analytical,X)
        f = -mu*div(grad(U)) - (mu+lamda)*grad(div(U)) - grad(P)
	L = -inner(f,v)*dx#-inner(v,p_analytical*n)*ds(2)# - inner(v,p_analytical*n)*ds(2) -inner(v,p_analytical*n)*ds(3) - inner(v,p_analytical*n)*ds(4) - inner(f,v)*dx# + k*inner(q,dot(grad(P),n))*ds(3)# + k*inner(q,dot(grad(P),n))*ds(2) + k*inner(q,dot(grad(P),n))*ds(3) + k*inner(q,dot(grad(P),n))*ds(4)

	A = a+L

	# Set up boundary condition (apply strong BCs)
	noslip = Constant((0.0,0.0))


	bcuT = DirichletBC(VQ.sub(0),u_analytical,top)
	bcuB = DirichletBC(VQ.sub(0),u_analytical,bottom)
        bcu_in = DirichletBC(VQ.sub(0),u_analytical,left)
        bcp = DirichletBC(VQ.sub(1),p_analytical,top)
        bcp2 = DirichletBC(VQ.sub(1),p_analytical,right)
        bcu_out = DirichletBC(VQ.sub(0),u_analytical,right)
        bcp3 = DirichletBC(VQ.sub(1),p_analytical,bottom)
        bcp4 = DirichletBC(VQ.sub(1),p_analytical,left)
	bcs = [bcuT,bcuB,bcu_in,bcu_out,bcp,bcp2,bcp3,bcp4]

	# Solution function
	UP_ = Function(VQ)

	# Assemble and apply boundary conditions
	"""
	A = assemble(a)
	b = assemble(L)
	[bc.apply(A, b) for bc in bcs]
	"""


	# Solve system


	solve(lhs(A) == rhs(A), UP_, bcs=bcs)
	u,p = UP_.split(True)
	e1 = errornorm(u,interpolate(u_analytical,W))
	ep1 = errornorm(p,interpolate(p_analytical,X))
	h1 = mesh.hmin()
	print 'N: ',N,' error u: %g ' %e1, 'rate : %g' %(ln(e1/e0)/ln(h1/h0))
	print 'N: ',N,' error p: %g ' %ep1, 'rate : %g' %(ln(ep1/ep0)/ln(h1/h0))
	#plot(p)
	e0 = e1
	ep0 = ep1
	h0 = mesh.hmin()
	#interactive()

ufile = File('res/u.pvd')
pfile = File('res/p.pvd')
ufile << u
pfile << p

plot(u)
plot(p)
plot(U)
plot(P)
interactive()
