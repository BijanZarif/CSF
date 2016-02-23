from domain import *
from pylab import zeros, where
import pylab as plt
import sys
set_log_active(False)
parameters['allow_extrapolation']=True

ufile = File("results_mono/velocity.pvd") # xdmf
pfile = File("results_mono/pressure.pvd")
dfile = File("results_mono/dU.pvd")
tfile = File("results_mono/U.pvd")


initial_files = False
if initial_files:
	t_folder = 20009998
	folder = 'initial_data_t_%d' % t_folder
	mesh = Mesh('%s/mesh.xml'%folder)
	SD = MeshFunction('size_t',mesh,'%s/mesh_func.xml'%folder)
	boundaries = MeshFunction('size_t',mesh,'%s/facet_func.xml'%folder)	
	V = VectorFunctionSpace(mesh,'CG',2)
	u1 = Function(V,'%s/u.xml'%folder)
	eta = Function(V,'%s/U.xml'%folder)
	d1 = Function(V, '%s/d.xml'%folder)
else:
	mesh = UnitSquareMesh(20,10)
	V = VectorFunctionSpace(mesh,'CG',2)

	SD = MeshFunction('uint', mesh, mesh.topology().dim())
	SD.set_all(0)
	#Solid().mark(SD,1)

	# DEFINING BOUNDARIES
	boundaries = FacetFunction("size_t",mesh)
	boundaries.set_all(0)

	class inlet(SubDomain):
		def inside(self,x,on_bnd):
			return on_bnd and near(x[0],0)

	class Top(SubDomain):
		def inside(self,x,on_bnd):
			return on_bnd and near(x[1],1)
	
	class Bottom(SubDomain):
		def inside(self,x,on_bnd):
			return on_bnd and near(x[1],0)

	class outlet(SubDomain):
		def inside(self,x,on_bnd):
			return on_bnd and near(x[0],1.0)

	inlet().mark(boundaries,1)
	outlet().mark(boundaries,2)
	Top().mark(boundaries,3)	
	Bottom().mark(boundaries,4)
	u1 = Function(V)#,'initial_data/u.xml')
	eta = Function(V)
	d1 = Function(V)

u0 = Function(V)


dt = 0.05   # use 0.0003 for oscillations
T = 10
# TEST AND TRIALFUNCTIONS

P = FunctionSpace(mesh,'CG',1)
W = FunctionSpace(mesh,'CG', 2)
VPW = MixedFunctionSpace([V,P,V])

v,p,w = TrialFunctions(VPW)
phi,eta,psi = TestFunctions(VPW)

# PHYSICAL PARAMETERS
h = mesh.hmin()

# FLUID
FSI = 1

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
noslip = Constant((0.0,0.0))
vf = Expression(('0.5*(1-cos(t*pi/2))*1.5*U*x[1]*(H-x[1])/pow((H/2.0),2)','0.0'),H=1,U=U_,t=2)

slipx = Expression(('speed','0'),speed=dt*0.1)
slipy = Expression(('0','-speed*x[0]'),speed=dt*0.1)


bcv1 = DirichletBC(VPW.sub(0),vf,boundaries,1)     # inlet
bcv3 = DirichletBC(VPW.sub(0),slipy,boundaries,3)
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4)
bcv = [bcv1,bcv3,bcv4]

# SOLID

# MESH DISPLACEMENT

bcw1 = DirichletBC(VPW.sub(2),noslip,boundaries,1)  # inlet
bcw2 = DirichletBC(VPW.sub(2),noslip,boundaries,2)
bcw3 = DirichletBC(VPW.sub(2),slipy,boundaries,3)
bcw4 = DirichletBC(VPW.sub(2),noslip,boundaries,4)
bcw = [bcw3,bcw4]


# CREATE FUNCTIONS

v1 = Function(V)#,'initial_data/u.xml')
U = Function(V)
v0 = Function(V)


VPW_ = Function(VPW)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))
n = FacetNormal(mesh)
g = Constant((0,2.0))


dS = Measure('dS')[boundaries]
dx = Measure('dx')[SD]
ds = Measure('ds')[boundaries]

#print assemble(n[0]*ds(6)), 'negative if n out of circle'
#print assemble(n('-')[0]*dS(5)), 'positive if n out of solid'
#sys.exit()
# n('+') out of solid, n('-') into solid
# n into circle


def sigma_dev(U):
	return 2*mu_s*sym(grad(U)) + lamda*tr(sym(grad(U)))*Identity(2)

def sigma_f(v):
	return 2*mu_f*sym(grad(v))

epsilon = 1e8

aMS = rho_s/k*inner(v,phi)*dx(1,subdomain_data=SD) + \
	k*inner(sigma_dev(v),grad(phi))*dx(1,subdomain_data=SD)

LMS = rho_s/k*inner(v1,phi)*dx(1,subdomain_data=SD) - \
	inner(sigma_dev(U),grad(phi))*dx(1,subdomain_data=SD)# - \
	#inner(g,phi)*dx(1,subdomain_data=SD)

#AV = epsilon*inner(grad(u),grad(w))*dx(1,subdomain_data=SD) - epsilon*inner(grad(d),grad(w))*dx(1,subdomain_data=SD)
aDS = epsilon*inner(v,psi)*dx(1,subdomain_data=SD) - epsilon*inner(w,psi)*dx(1,subdomain_data=SD)


aS = aMS + aDS
LS = LMS

# FLUID

aMF =rho_f/k*inner(v,phi)*dx(0,subdomain_data=SD)+ \
	rho_f*inner(grad(v0)*(v-w),phi)*dx(0,subdomain_data=SD) - \
	 inner(p,div(phi))*dx(0,subdomain_data=SD) + \
	mu_f*inner(grad(v),grad(phi))*dx(0,subdomain_data=SD)

LMF = rho_f/k*inner(v1,phi)*dx(0,subdomain_data=SD)
'''
aMF = rho_f/k*inner(v,phi)*dx(0,subdomain_data=SD) + \
	rho_f*inner(grad(v0)*(v-w),phi)*dx(0,subdomain_data=SD) - \
	 inner(p,div(phi))*dx(0,subdomain_data=SD) + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx(0,subdomain_data=SD)

LMF = rho_f/k*inner(v1,phi)*dx(0,subdomain_data=SD)
'''
aDF = k*inner(grad(w),grad(psi))*dx(0,subdomain_data=SD)
LDF = -inner(grad(U),grad(psi))*dx(0,subdomain_data=SD)

aCF = -inner(div(v),eta)*dx(0,subdomain_data=SD)
LFQ = -Constant(0)*eta*dx(0,subdomain_data=SD)

aF = aMF + aDF + aCF
LF = LMF + LDF + LFQ

a = aS+aF
L = LS+LF
if initial_files:
	t = t_folder*1e-6
else:
	t = dt
count = 0




while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):
	b = assemble(L)
	if t > 10.0:
		bcw2.speed = k*0.1
	eps = 10
	k_iter = 0
	max_iter = 7
	while eps > 1E-6 and k_iter < max_iter:
	    A = assemble(a)
	    A.ident_zeros()
	    [bc.apply(A,b) for bc in bcv]
	    [bc.apply(A,b) for bc in bcw]
	    solve(A,VPW_.vector(),b,'lu')
	    v_,p_,w_ = VPW_.split(True)
	    eps = errornorm(v_,v0,degree_rise=3)
	    k_iter += 1
	    print 'k: ',k_iter, 'error: %.3e' %eps
	    v0.assign(v_)
	
	if count%1==0:
		ufile << v_
		pfile << p_
		dfile << w_
		tfile << U
	#w_.vector()[:] *= float(k)
	U.vector()[:] += float(k)*w_.vector()[:]
	mesh.move(w_)
	mesh.bounding_box_tree().build(mesh)
	# Move to next time step
	v1.assign(v_)
        print 't=%.4f'%t
	t += dt
	#print '%10d %10.4e %5.4e' %(LEN,xA-mesh.coordinates()[coord,0], yA-mesh.coordinates()[coord,1])
	count += 1


folder = 'initial_data_t_%d'% int(t*10**6)
File('%s/u.xml'%folder) << v_
File('%s/p.xml'%folder) << p_
File('%s/d.xml'%folder) << w_
File('%s/U.xml'%folder) << U
File('%s/mesh.xml'%folder) << mesh
File('%s/mesh_func.xml'%folder) << SD
File('%s/facet_func.xml'%folder) << boundaries   

