from domains_and_boundaries import *
set_log_active(False)
parameters['allow_extrapolation']=True

# REFINEMENT LEVEL AND YOUNGS MODULUS
RL = '1'
EPa = '62500'			# [Pa]
ref_level = int(RL)
E = Constant(float(EPa))

vfile = File("RESULTS/ref%s_E%s/velocity.pvd"%(RL,EPa)) # xdmf
pfile = File("RESULTS/ref%s_E%s/pressure.pvd"%(RL,EPa))
wfile = File("RESULTS/ref%s_E%s/dU.pvd"%(RL,EPa))
Ufile = File("RESULTS/ref%s_E%s/U.pvd"%(RL,EPa))

# PHYSICAL PARAMETERS
rho_f = Constant(1./1000)		# [g/mm]
nu_f = Constant(0.658)			# [mm**2/s]
mu_f = Constant(nu_f*rho_f)		# [g/(mm*s)]

Pr = 0.479
rho_s = Constant(1.75*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))

# MESH SETUP
Nx = 18*ref_level
Ny = 30*ref_level
P1 = Point(-9,0)
P2 = Point(9,60)
mesh = RectangleMesh(P1,P2,Nx,Ny)

SD = MeshFunction('size_t', mesh, mesh.topology().dim())
SD.set_all(1)
Fluid().mark(SD,0)
CSC().mark(SD,0)

# DEFINING BOUNDARIES
boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
Fluid_in().mark(boundaries,1)
Solid_in().mark(boundaries,2)
Fluid_out().mark(boundaries,3)
Solid_out().mark(boundaries,4)
Interface().mark(boundaries,5)
Fluid_walls().mark(boundaries,6)
CSC_bnd().mark(boundaries,7)
Fluid_in().mark(boundaries,1)
Fluid_out().mark(boundaries,3)

# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = VectorFunctionSpace(mesh,'CG', 1)
VPW = MixedFunctionSpace([V,P,W])

v,p,w = TrialFunctions(VPW)
phi,eta,psi = TestFunctions(VPW)

# PARAMETERS
dt = 0.002
T = 10
hmesh = mesh.hmin()
k = Constant(dt)
delta = 1e-8

# INITIAL AND BOUNDARY CONDITIONS
# FLUID
noslip = Constant((0.0,0.0))

pressure = Expression(('amp*sin(2*pi*t)'),t=0,amp=1)

bcv2 = DirichletBC(VPW.sub(0),noslip,boundaries,2) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Solid out
bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Fluid walls

bcv = [bcv2, bcv4, bcv6]

# SOLID
bcw1 = DirichletBC(VPW.sub(2),noslip,boundaries,1)  # Fluid in
bcw2 = DirichletBC(VPW.sub(2),noslip,boundaries,2)  # Solid in
bcw3 = DirichletBC(VPW.sub(2),noslip,boundaries,3)  # Fluid out
bcw4 = DirichletBC(VPW.sub(2),noslip,boundaries,4)  # Solid out
bcw6 = DirichletBC(VPW.sub(2),noslip,boundaries,6) # Fluid walls
bcw = [bcw1,bcw2,bcw3,bcw4,bcw6]

bcs = bcv+bcw

# CREATE FUNCTIONS
v0 = Function(V)
v1 = Function(V)
U = Function(W)

VPW_ = Function(VPW)


# MEASUREMENTS FOR INTEGRATION
dS = Measure('dS')[boundaries]
dx = Measure('dx')[SD]
ds = Measure('ds')[boundaries]

dx_f = dx(0,subdomain_data=SD)
dx_s = dx(1,subdomain_data=SD)

k = Constant(dt)
n = FacetNormal(mesh)

# CONSTITUTIVE RELATIONSHIP
def sigma_dev(U):
	return 2*mu_s*sym(grad(U)) + lamda*tr(sym(grad(U)))*Identity(2)

# SET UP EQUATIONS
# SOLID
aMS = rho_s/k*inner(v,phi)*dx_s \
	+ rho_s*inner(grad(v0)*v,phi)*dx_s \
	+ inner(sigma_dev(k*v),grad(phi))*dx_s

LMS = rho_s/k*inner(v1,phi)*dx_s \
	+ inner(sigma_dev(U),grad(phi))*dx_s

aDS = 1/delta*inner(v,psi)*dx_s - 1/delta*inner(w,psi)*dx_s

aS = aMS + aDS
LS = LMS

# FLUID
penalty = 0.01*hmesh
aMF = rho_f/k*inner(v,phi)*dx_f \
	+ rho_f*inner(grad(v0)*(v-w),phi)*dx_f \
	- inner(p,div(phi))*dx_f \
	+ 2*mu_f*inner(sym(grad(v)),grad(phi))*dx_f \
	- mu_f*inner(grad(v).T*n,phi)*ds(1) \
	- mu_f*inner(grad(v).T*n,phi)*ds(3) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(1) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(3)

LMF = rho_f/k*inner(v1,phi)*dx_f - \
	inner(pressure*n,phi)*ds(1) + \
	inner(pressure*n,phi)*ds(3)

aDF = k*inner(grad(w),grad(psi))*dx_f \
	+ k*inner(grad(w('-'))*n('-'),psi('-'))*dS(5)
LDF = -inner(grad(U),grad(psi))*dx_f \
	+ inner(grad(U('-'))*n('-'),psi('-'))*dS(5)

aCF = -inner(div(v),eta)*dx_f

aF = aMF + aDF + aCF
LF = LMF + LDF

# LINEAR AND BILINEAR FORMS
a = aS+aF
L = LS+LF

t = dt
count = 0
while t < T + DOLFIN_EPS:
	pressure.t = t
	b = assemble(L)
	eps = 10			# parameters for Picard iteration
	k_iter = 0
	max_iter = 5
	while eps > 1E-6 and k_iter < max_iter:
		A = assemble(a)
		A.ident_zeros()
		[bc.apply(A,b) for bc in bcs]
		solve(A,VPW_.vector(),b,'lu')
		v_,p_,w_ = VPW_.split(True)
		eps = errornorm(v_,v0,degree_rise=3)
		k_iter += 1
		print 'k: ',k_iter, 'error: %.3e' %eps
		v0.assign(v_)
	if count%5==0: # Save every fifth state
		vfile << v_
		pfile << p_
		wfile << w_
		Ufile << U

	# Update mesh and move to next time step
	w_.vector()[:] *= float(k)
	U.vector()[:] += w_.vector()[:]
	ALE.move(mesh,w_)
	mesh.bounding_box_tree().build(mesh)
	v1.assign(v_)
	print 't=%.4f'%t
	t += dt
	count += 1
