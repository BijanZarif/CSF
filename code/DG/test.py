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

class Interface(SubDomain):
        def inside(self, x, on_boundary):
                return abs(x[1] - 0.3) < DOLFIN_EPS and abs(x[0]-0.4) < DOLFIN_EPS and abs(x[0]-0.6) < DOLFIN_EPS

class Domain2(SubDomain):
        def inside(self, x, on_boundary):
                return (x[1]<0.3+DOLFIN_EPS) and abs(x[0]-0.5) < 0.2

h0 = 1
e0 = 1
ep0 = 1

for N in [20]:#,8,16,32,64]:
	mesh = UnitSquareMesh(N,N)
	right = Right()
	left = Left()
	top = Top()
	bottom = Bottom()
        interface = Interface()
        domain2 = Domain2()
        
        subdomains = MeshFunction('uint',mesh,mesh.topology().dim())
        subdomains.set_all(0)
        domain2.mark(subdomains,1)
	mesh = UnitSquareMesh(N,N)
	boundaries = FacetFunction('size_t',mesh)
	boundaries.set_all(0)
	left.mark(boundaries,1)
	right.mark(boundaries,2)
	top.mark(boundaries,3)
	bottom.mark(boundaries,3)

	ds = Measure('ds')[boundaries]
        dx = Measure('dx')[subdomains]

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

	# Source term
	f = Constant((0.0,0.0))

	# Penalty term
	alpha = Constant(14.0)
	beta = Constant(10.0)

	# Mesh-related functions
	n = FacetNormal(mesh)
	h = CellSize(mesh)
	h_avg = (h('+') + h('-'))/2

	# ( dot(v, n) + |dot(v, n)| )/2.0
	p_in = Constant(-2.0)
	p_out = Constant(0.0)

	n = FacetNormal(mesh)

	def eps(u):
		return 0.5*(grad(u) + grad(u).T)

	# Bilinear form
	a_int = 2*mu*inner(eps(u),eps(v))*dx + inner(p,div(v))*dx + inner(q,div(u))*dx


	a2_int = mu*inner(grad(u),grad(v))*dx(0,subdomain_data=subdomains) + inner(p,div(v))*dx(0, subdomain_data=subdomains)+ inner(q,div(u))*dx(0, subdomain_data=subdomains) - 1.0/N**2*inner(grad(p),grad(q))*dx(0,subdomain_data=subdomains)
        a2_int2 = mu*inner(grad(u),grad(v))*dx(1,subdomain_data=subdomains) + inner(p,div(v))*dx(1, subdomain_data=subdomains)+ inner(q,div(u))*dx(1, subdomain_data=subdomains) - 1.0/N**2*inner(grad(p),grad(q))*dx(1,subdomain_data=subdomains)
        

        
        
	a2_diff = - mu('+')*inner(dot(avg(grad(u)),n('+')),jump(v))*dS \
                  - mu('+')*inner(dot(avg(grad(v)),n('+')),jump(u))*dS \
		  + mu('+')*beta('+')/h('+')*inner(jump(u),jump(v))*dS \
		  - inner(avg(p),dot(jump(v),n('+')))*dS \
		  - inner(avg(q),dot(jump(u),n('+')))*dS \
                  + beta('+')/h('+')*inner(jump(q),jump(p))*dS \
		  - mu*inner(dot(grad(u),n),v)*ds(3) \
		  - mu*inner(dot(grad(v),n),u)*ds(3) \
		  + mu*alpha/h*inner(u,v)*ds(3)
        


	a_diff =  - 2*mu('+')*inner(dot(avg(eps(u)),n('+')),jump(v))*dS \
                  - 2*mu('+')*inner(dot(avg(eps(v)),n('+')),jump(u))*dS \
		  + 2*mu('+')*beta('+')/h('+')*inner(jump(u),jump(v))*dS \
		  - inner(avg(p),dot(jump(v),n('+')))*dS \
		  - inner(avg(q), dot(jump(u),n('+')))*dS \
		  - 2*mu*inner(dot(eps(u),n),v)*ds(3) \
		  - 2*mu*inner(dot(eps(v),n),u)*ds(3) \
		  + 2*mu*alpha/h*inner(u,v)*ds(3) \
		  - mu*inner(dot(grad(u).T,n),v)*ds(1) \
		  - mu*inner(dot(grad(u).T,n),v)*ds(2)


	#a = a_int + a_diff# + a_fac + a_vel
        a = a2_int + a2_int2 #a2_diff

	# Linear form
	L = -inner(v,p_in*n)*ds(1,subdomain_data = subdomains) - inner(v,p_out*n)*ds(2,subdomain_data = subdomains) + inner(f,v)*dx(1,subdomain_data=subdomains)

	A = a+L

	# Set up boundary condition (apply strong BCs)
	noslip = Constant((0.0,0.0))
	u_analytical = Expression(('x[1]*(1-x[1])','0.0'))
	p_analytical = Expression('-2+2*x[0]')

	bcuT = DirichletBC(VQ.sub(0),noslip,top)
	bcuB = DirichletBC(VQ.sub(0),noslip,bottom)
        bcu_in = DirichletBC(VQ.sub(0),u_analytical,left)
	bcs = [bcuT, bcuB]

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
interactive()
