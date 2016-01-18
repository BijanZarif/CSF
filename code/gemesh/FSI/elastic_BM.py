from domains_and_boundaries import *
from pylab import zeros, where
import pylab as plt
import sys
set_log_active(False)
parameters['allow_extrapolation']=True

ufile = File("results_mono/velocity.pvd")
pfile = File("results_mono/pressure.pvd")
dfile = File("results_mono/dU.pvd")
tfile = File("results_mono/U.pvd")

initial_files = False
if initial_files:
    t_folder = 1500000
    folder = 'initial_data_t_%d' % t_folder

dt = 0.001
T = 20*dt
mesh = Mesh('../meshes/cord_w_csc_mm.xml')

#mesh = refine(mesh)
#mesh = refine(mesh)
LEN = len(mesh.coordinates())
print 'len(mesh.coord)', LEN
SD = MeshFunction('uint', mesh, mesh.topology().dim())

SD.set_all(1)
cord().mark(SD,0)
inside_csc().mark(SD,0)

# DEFINING BOUNDARIES
boundaries = FacetFunction("size_t",mesh)
outer_wall().mark(boundaries,1)
inner_wall().mark(boundaries,2)
top_outer().mark(boundaries,3)
bottom_outer().mark(boundaries,4)
bottom_inner().mark(boundaries,5)
top_inner().mark(boundaries,6)
csc_boundary().mark(boundaries,7)

# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)
W = FunctionSpace(mesh,'CG', 2)
VQW = MixedFunctionSpace([V,Q,V])

u,p,d = TrialFunctions(VQW)
v,q,w = TestFunctions(VQW)

# PHYSICAL PARAMETERS
h = mesh.hmin()

# FLUID
FSI = 2

nu = Constant(0.658)
rho_f = Constant(1.0*1e-6)
mu_f = Constant(rho_f*nu)

U_ = [0.2, 1.0, 2.0][FSI-1]   # Reynolds vel
# SOLID
E = 5.0*10**3
Pr = 0.479

rho_s = 1.0
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))

K_Perm = 1.4*10**-6
kappa = K_perm/mu_f


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0.0,0.0))
vf = Expression(('0,0','1.5*U*x[1]*(H-x[1])/pow((H/2.0),2)'),H=0.075,U=U_)

p_in = Expression('100*t*cos(4*pi*t)',t=0)
p_out = Expression('cos(4*pi*t)',t=0)

bcu1 = DirichletBC(VQW.sub(0),noslip,boundaries,1) 	# walls
bcu2 = DirichletBC(VQW.sub(0),noslip,boundaries,2) 	# interface
bcu3 = DirichletBC(VQW.sub(0),vf,boundaries,3)		# fluid inlet
bcu4 = DirichletBC(VQW.sub(0),noslip,boundaries,4) 	# fluid outlet
bcu5 = DirichletBC(VQW.sub(0),noslip,boundaries,5) 	# elastic outlet
bcu6 = DirichletBC(VQW.sub(0),noslip,boundaries,6)      # elastic inlet
bcu7 = DirichletBC(VQW.sub(0),noslip,boundaries,7)      # csc boundary

bcu = [bcu1,bcu6,bcu5]

bcp3 = DirichletBC(VQW.sub(1),p_in,boundaries,3)
bcp6 = DirichletBC(VQW.sub(1),p_in,boundaries,6)
bcp4 = DirichletBC(VQW.sub(1),p_out,boundaries,4)
bcp5 = DirichletBC(VQW.sub(1),p_out,boundaries,5)

bcp = [bcp6,bcp5,bcp3,bcp4]

# SOLID

# MESH DISPLACEMENT

bcU1 = DirichletBC(VQW.sub(2),noslip,boundaries,1)  # walls
bcU2 = DirichletBC(VQW.sub(2),noslip,boundaries,2)  # interface
bcU3 = DirichletBC(VQW.sub(2),noslip,boundaries,3)  # fluid inlet
bcU4 = DirichletBC(VQW.sub(2),noslip,boundaries,4)  # fluid outlet
bcU5 = DirichletBC(VQW.sub(2),noslip,boundaries,5)  # elastic outlet
bcU6 = DirichletBC(VQW.sub(2),noslip,boundaries,6)  # elastic inlet

bcP = [bcU1]#,bcU3,bcU4,bcU5,bcU6]

# CREATE FUNCTIONS
u0 = Function(V)


UPR1 = Function(VQW)
u1,p1,d1 = UPR1.split(True)
UPR = Function(VQW)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))
n = FacetNormal(mesh)
g = Constant((0,2.0))


dS = Measure('dS')[boundaries]
dx = Measure('dx')[SD]
ds = Measure('ds')[boundaries]

#print assemble(n('-')[0]*dS(5))
# n('-') into solid, n('+') out of solid
# n into circle

def sigma_f(u,p):
    return -p*Identity(2) + 2*mu_f*sym(grad(u))

def sigma_dev(d):
    return 2*mu_s*sym(grad(d)) + lamda*tr(sym(grad(d)))*Identity(2)

# FLUID

aF = rho_f/k*inner(u,v)*dx(1,subdomain_data=SD) + \
     rho_f*inner(grad(u0)*u,v)*dx(1,subdomain_data=SD) + \
     mu_f*inner(grad(u),grad(v))*dx(1,subdomain_data=SD) - \
     inner(p,div(v))*dx(1,subdomain_data=SD) - \
     inner(q,div(u))*dx(1,subdomain_data=SD) + \
     inner(p('-')*n('-'),v('-'))*dS(2)

LF = rho_f/k*inner(u1,v)*dx(1,subdomain_data=SD) - \
     inner(p_in*n,v)*ds(3) + \
     inner(sigma_dev(d1('-'))*n('-'),v('-'))*dS(2)

# SOLID
aS = 1./k*inner(u,w)*dx(0,subdomain_data=SD) - \
     inner(p,div(w))*dx(0,subdomain_data=SD) - \
     inner(w('+'),sigma_f(u('+'),p('+'))*n('+'))*dS(2)

LS = 1./k*inner(u1,w)*dx(0,subdomain_data=SD) - \
     inner(grad(w),sigma_dev(d1))*dx(0,subdomain_data=SD)

aU = 1./k*inner(d,v)*dx(0,subdomain_data=SD) - inner(u,v)*dx(0,subdomain_data=SD)
LU = 1./k*inner(d1,v)*dx(0,subdomain_data=SD)


aV = inner(u,w)*dx(0,subdomain_data=SD) - 1./k*inner(d,w)*dx(0,subdomain_data=SD) + kappa/mu_f*inner(grad(p),w)*dx(0,subdomain_data=SD)
LV = -1./k*inner(d1,w)*dx(0,subdomain_data=SD)


#POISSON
aP = inner(grad(d),grad(w))*dx(1,subdomain_data=SD) 
LP = inner(f,w)*dx(1,subdomain_data=SD)


a = aF + aS + aP + aU# + aV
L = LF + LS + LP + LU# + LV

if initial_files:
	t = t_folder*1e-6

else:
	t = dt

change_vf = True


xA = where(mesh.coordinates()[:,0] == 0.6)
yA = where(mesh.coordinates()[:,1] == 0.2)

n1 = as_vector((1.0,0))
n2 = as_vector((0,1.0))
nx = dot(n1,n)
ny = dot(n2,n)
nt = as_vector((ny,-nx))
xA = plt.where(abs(mesh.coordinates()[:,0] - 0.6) < h)
yA = plt.where(abs(mesh.coordinates()[:,1] - 0.2) < h)
for idx in xA[0]:
    if idx in yA[0]:
        coord = idx
xA = float(mesh.coordinates()[coord,0])
yA = float(mesh.coordinates()[coord,1])
count = 0

while t < T + DOLFIN_EPS:
        u1,p1,d1 = UPR1.split(True)
	p_in.t = t
	#p_out.t = t
	
	b = assemble(L)
	eps = 10
	k_iter = 0
	max_iter = 5
	while eps > 1E-4 and k_iter < max_iter:
	    A = assemble(a)
	    A.ident_zeros()
	    [bc.apply(A,b) for bc in bcu]
	    [bc.apply(A,b) for bc in bcP]
	    solve(A,UPR.vector(),b,'lu')
	    u_,p_,d_ = UPR.split(True)
	    eps = errornorm(u_,u0,degree_rise=3)
	    k_iter += 1
	    print 'k: ',k_iter, 'error: %.3e' %eps
	    u0.assign(u_)


	#bcu[-1] = DirichletBC(VQW.sub(0),d_,boundaries,5) # interface
        #D = interpolate(d_-d1,V)
	ufile << u_
	pfile << p_
	dfile << d_
        d1.vector()[:] += u_.vector()[:]
	mesh.move(d_)
	mesh.bounding_box_tree().build(mesh)
	# Move to next time step
        UPR1.assign(UPR)
        u1,p1,d1 = UPR1.split()
        print 't=%.4f'%t
	t += dt

	count += 1





folder = 'initial_data_t_%d'% int(t*10**6)
File('%s/u.xml'%folder) << u_
File('%s/p.xml'%folder) << p_
File('%s/d.xml'%folder) << d_
    
