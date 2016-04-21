from domains_and_boundaries import *
from numpy import zeros, where
import sys
#set_log_active(False)
parameters['allow_extrapolation']=True

ufile = File("results_mono/velocity.pvd") # xdmf
pfile = File("results_mono/pressure.pvd")
dfile = File("results_mono/dU.pvd")
tfile = File("results_mono/U.pvd")

ref_level = 1
## from Karen/Nina
E = (5*10**3) # Pa --- Ozawa rabbit spinal cord 5*10**3 (No Pia) 16*10^3 (w Pia)
Pr = 0.479#0.479       #


rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)


## from Karen/Nina
E = (5*10**3) # Pa --- Ozawa rabbit spinal cord 5*10**3 (No Pia) 16*10^3 (w Pia)
Pr = 0.479#0.479       #

rho_s = Constant(1.75*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))


Nx = 18*ref_level
Ny = 30*ref_level

P1 = Point(-9,0)
P2 = Point(9,60)
mesh = RectangleMesh(P1,P2,18,30)
LEN = len(mesh.coordinates())
print 'len(mesh.coord)', LEN

SD = MeshFunction('size_t', mesh, mesh.topology().dim())
SD.set_all(0)
Solid().mark(SD,1)
CSC().mark(SD,1)

# DEFINING BOUNDARIES
boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
Fluid_in_l().mark(boundaries,1)
Fluid_in_r().mark(boundaries,2)
Solid_in().mark(boundaries,3)
Fluid_out().mark(boundaries,4)
Solid_out().mark(boundaries,5)
Interface().mark(boundaries,6)
Fluid_walls().mark(boundaries,7)
CSC_bnd().mark(boundaries,8)


dt = 0.002   # use 0.0003 for oscillations
T = 10
# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = VectorFunctionSpace(mesh,'CG', 1)
VPW = MixedFunctionSpace([V,P,W])

vpw = TrialFunction(VPW)
PhiEtaPsi = TestFunction(VPW)

vpw_ = Function(VPW)
vpw1 = Function(VPW)

v,p,w = split(vpw)
v1,p1,w1 = split(vpw1)
phi,eta,psi = split(PhiEtaPsi)

# PHYSICAL PARAMETERS
h = mesh.hmin()

# FLUID
FSI = 3

nu = 10**-3
rho_f = 1.0*1e3
mu_f = rho_f*nu
U_ = [0.2, 1.0, 2.0][FSI-1]   # Reynolds vel

# SOLID
Pr = 0.4
mu_s = [0.5, 0.5, 2.0][FSI-1]*1e6
rho_s = [1.0, 10, 1.0][FSI-1]*1e3
lamda = 2*mu_s*Pr/(1-2.*Pr)


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0,0))

bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Fluid out
bcv5 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Solid out
bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Interface
bcv7 = DirichletBC(VPW.sub(0),noslip,boundaries,7) # Fluid walls


bcv = [bcv3, bcv5, bcv7] # don't use bcv6 for FSI

pressure = Expression(('amp*sin(2*pi*t)'),t=0,amp=1)

# SOLID

# MESH DISPLACEMENT

bcw1 = DirichletBC(VPW.sub(2),noslip,boundaries,1)  # Fluid in_l
bcw2 = DirichletBC(VPW.sub(2),noslip,boundaries,2)  # Fluid in_r
bcw3 = DirichletBC(VPW.sub(2),noslip,boundaries,3)  # Solid in
bcw4 = DirichletBC(VPW.sub(2),noslip,boundaries,4)  # Fluid out
bcw5 = DirichletBC(VPW.sub(2),noslip,boundaries,5)  # Solid out
bcw6 = DirichletBC(VPW.sub(2),noslip,boundaries,6)  # Interface
bcw7 = DirichletBC(VPW.sub(2),noslip,boundaries,7) # Fluid walls
bcw = [bcw1,bcw2,bcw3,bcw4,bcw5,bcw7]

# CREATE FUNCTIONS
	
U = Function(W)
U0 = Function(W)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))
n = FacetNormal(mesh)
g = Constant((0,2.0))


dS = Measure('dS')[boundaries]
dx = Measure('dx')[SD]
ds = Measure('ds')[boundaries]


def Eij(U,U0):
	return sym(grad(U)) + dot(grad(U),grad(U0))

def sigma_dev(U):
	return 2*mu_s*sym(grad(U)) + lamda*tr(sym(grad(U)))*Identity(2)

def sigma_2(U,U0):
	return 2*mu_s*Eij(U,U0) + lamda*tr(Eij(U,U0))*Identity(2)


epsilon = 1e8

U_ = U + k*v      # discretize in time here

dx_s = dx(1,subdomain_data=SD)
dx_f = dx(0,subdomain_data=SD)

D = inner(sigma_dev(U_),grad(phi))*dx_s

aMS = rho_s/k*inner(v,phi)*dx_s \
	+ rho_s*inner(grad(v)*v,phi)*dx_s \
	+ lhs(D)

LMS = rho_s/k*inner(v1,phi)*dx_s \
	+ rhs(D)

aDS = epsilon*inner(v,psi)*dx_s - epsilon*inner(w,psi)*dx_s

aCS = epsilon**-1*inner(p,eta)*dx_s

aS = aMS + aDS + aCS
LS = LMS 

penalty = 0.05*mesh.hmin()
# FLUID

aMF = rho_f/k*inner(v,phi)*dx_f \
	+ rho_f*inner(grad(v)*(v-w),phi)*dx_f \
	+ inner(p,div(phi))*dx_f \
	+ 2*mu_f*inner(sym(grad(v)),grad(phi))*dx_f \
	- mu_f*inner(grad(v).T*n,phi)*ds(1) \
	- mu_f*inner(grad(v).T*n,phi)*ds(2) \
	- mu_f*inner(grad(v).T*n,phi)*ds(4) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(1) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(2) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(4)


LMF = rho_f/k*inner(v1,phi)*dx_f - \
	inner(pressure*n,phi)*ds(1) - \
	inner(pressure*n,phi)*ds(2) + \
	inner(pressure*n,phi)*ds(4)

aDF = k*inner(grad(w),grad(psi))*dx_f
LDF = -inner(grad(U),grad(psi))*dx_f

aFQ = -inner(div(v),eta)*dx_f


aF = aMF + aDF + aFQ
LF = LMF + LDF

a = aS+aF
L = LS+LF

F = a-L
F = action(F,vpw_)
bcs = bcv+bcw

count = 0
t = dt

while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):
	if t <= 2.0:
		pressure.amp = 5*t
	else:
		pressure.amp = 10
	pressure.t = t

	solve(F==0,vpw_,bcs,solver_parameters={"newton_solver": {"relative_tolerance": 1e-6}})
	v_,p_,w_ = vpw_.split(True)
	U0.assign(U)
	U0.vector()[:] += float(dt)*(w_.vector()[:])
	if count%20==0:
		ufile << v_
		pfile << p_
		dfile << w_
		tfile << U
	
	w_.vector()[:] *= float(k)
	U.vector()[:] += w_.vector()[:]
	mesh.move(w_)
	mesh.bounding_box_tree().build(mesh)
	vpw1.assign(vpw_)
	print 't=%.4f'%t
	t += dt

	count += 1


