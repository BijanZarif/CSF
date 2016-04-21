from dolfin import *
import sys
set_log_active(False)
ufile = File('res/u.pvd')
pfile = File('res/p.pvd')
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
'''
a2_int = mu*inner(grad(u),grad(v))*dx + inner(p,div(v))*dx + inner(q,div(u))*dx# - 1.0/N**2*inner(grad(p),grad(q))*dx

a2_diff = - mu('+')*inner(dot(avg(grad(u)),n('+')),jump(v))*dS \
              - mu('+')*inner(dot(avg(grad(v)),n('+')),jump(u))*dS \
		      + mu('+')*beta('+')/h('+')*inner(jump(u),jump(v))*dS \
		  	  - inner(avg(p),dot(jump(v),n('+')))*dS \
		   	  - inner(avg(q),dot(jump(u),n('+')))*dS \
              + beta('+')/h('+')*inner(jump(q),jump(p))*dS \
		  	  + mu*alpha/h*inner(u,v)*ds(3)
'''
for N in [16]:#[4,8,16,32,64]:
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
	V = VectorFunctionSpace(mesh, "CG",2)
	VQ = V*Q
	W = VectorFunctionSpace(mesh,'CG',3)
	X = FunctionSpace(mesh,'CG',3)
	# Create velocity Function

	print VQ.dim()

	# Test and trial functions
	v,q  = TestFunctions(VQ)
	u,p = TrialFunctions(VQ)
	# Diffusivity
	mu = Constant(1.0)
	rho = Constant(1.0)

	# Source term
	f = Constant((0.0,0.0))

	# Penalty term
	alpha = Constant(10)
	beta = Constant(10)

	# Mesh-related functions
	n = FacetNormal(mesh)
	h = CellSize(mesh)
	h_avg = (h('+') + h('-'))/2

	# ( dot(v, n) + |dot(v, n)| )/2.0
	p_in = Expression('2',t=0)
	p_out = Constant(0.0)
	dt = Constant(0.01)

	n = FacetNormal(mesh)

	def eps(u):
		return 0.5*(grad(u) + grad(u).T)
	
	u0 = Function(V)
	u1 = Function(V)
	# Bilinear form
	a_acc = rho/dt*inner(u,v)*dx + rho*inner(grad(u0)*u,v)*dx
	a_int = 2*mu*inner(eps(u),eps(v))*dx - inner(p,div(v))*dx - inner(q,div(u))*dx

	a_diff=-2*mu('+')*inner(avg(eps(u))*n('+'),jump(v))*dS \
          - 2*mu('+')*inner(avg(eps(v))*n('+'),jump(u))*dS \
		  + 2*mu('+')*beta('+')/h('+')*inner(jump(u),jump(v))*dS \
		  - inner(avg(p),dot(jump(v),n('+')))*dS \
		  - inner(avg(q),dot(jump(u),n('+')))*dS \
		  + 2*mu*alpha/h*inner(u,v)*ds(3) \
		  - mu*inner(grad(u).T*n,v)*ds(1) \
		  - mu*inner(grad(u).T*n,v)*ds(2)


	a = a_int + a_diff# + a_acc
	#a = a2_int + a2_diff

	# Linear form
	L = - inner(v,p_in*n)*ds(1) - inner(v,p_out*n)*ds(2) + rho/dt*inner(u1,v)*dx

	# Set up boundary condition (apply strong BCs)
	noslip = Constant((0.0,0.0))
	u_analytical = Expression(('x[1]*(1-x[1])','0.0'))
	p_analytical = Expression('2-2*x[0]')

	bcuT = DirichletBC(VQ.sub(0),noslip,boundaries,3)
	bcu_in = DirichletBC(VQ.sub(0),u_analytical,left)
	bcs = [bcuT]

	# Solution function
	UP_ = Function(VQ)

	# Assemble and apply boundary conditions
	
	

	T = 3
	t = float(dt)

	while T > t:
		p_in.t = t
		b = assemble(L)
		err = 10
		max_iter = 3
		k_iter = 0
		while err > 1e-2 or k_iter < max_iter:
			A = assemble(a)
			#[bc.apply(A, b) for bc in bcs]
	# Solve system
			solve(A, UP_.vector(),b)
			u_,p_ = UP_.split(True)
			err = errornorm(u_,u0)
			u0.assign(u_)
			k_iter += 1
	
		u1.assign(u_)
		e1 = errornorm(u_,interpolate(u_analytical,W))
		ep1 = errornorm(p_,interpolate(p_analytical,X))
		h1 = mesh.hmin()
		print 't=',t,'  error u = %2.4e   error p = %2.4e'%(e1,ep1)
		e0 = e1
		ep0 = ep1
		h0 = mesh.hmin()
		#interactive()
		t+=float(dt)
		ufile << u_
		pfile << p_

