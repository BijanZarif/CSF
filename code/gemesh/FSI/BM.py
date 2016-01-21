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
T = 10*dt
mesh = Mesh('../meshes/cord_w_csc_mm.xml')

#mesh = refine(mesh)
#mesh = refine(mesh)
LEN = len(mesh.coordinates())
print 'len(mesh.coord)', LEN
SD = MeshFunction('uint', mesh, mesh.topology().dim())

SD.set_all(0)
cord().mark(SD,1)
inside_csc().mark(SD,1)

# DEFINING BOUNDARIES
boundaries = FacetFunction("size_t",mesh)
outer_wall().mark(boundaries,1)
inner_wall().mark(boundaries,2)
top_right().mark(boundaries,3)
bottom_outer().mark(boundaries,4)
bottom_inner().mark(boundaries,5)
top_inner().mark(boundaries,6)
csc_boundary().mark(boundaries,7)
top_left().mark(boundaries,8)

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

nu = 0.658
rho_f = 1.0*1e-3
mu_f = rho_f*nu

# SOLID
E = 5*10**3
Pr = 0.479
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))

kappa = Constant(1.4*10**(-15)*(10**6))


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0.0,0.0))
u_in_R = Expression(('0,0','0.01*x[0]*(x[0]-x0)'),x0=x1)
u_in_L = Expression(('0,0','0.01*x[0]*(x0-x[0])'),x0=x1)

p_in = Expression('1000*cos(4*pi*t)',t=0)
p_out = Expression('1000*sin(4*pi*t)',t=0)

bcu1 = DirichletBC(VQW.sub(0),noslip,boundaries,1) 	# walls
bcu2 = DirichletBC(VQW.sub(0),noslip,boundaries,2) 	# interface
bcu3 = DirichletBC(VQW.sub(0),u_in_R,boundaries,3)		# fluid inlet right
bcu4 = DirichletBC(VQW.sub(0),noslip,boundaries,4) 	# fluid outlet
bcu5 = DirichletBC(VQW.sub(0),noslip,boundaries,5) 	# elastic outlet
bcu6 = DirichletBC(VQW.sub(0),noslip,boundaries,6)      # elastic inlet
bcu7 = DirichletBC(VQW.sub(0),noslip,boundaries,7)      # csc boundary
bcu8 = DirichletBC(VQW.sub(0),u_in_L,boundaries,8)    # fluid inlet Left

bcu = [bcu1,bcu3,bcu8,bcu6,bcu5]

bcp3 = DirichletBC(VQW.sub(1),p_in,boundaries,3)
bcp6 = DirichletBC(VQW.sub(1),p_in,boundaries,6)
bcp4 = DirichletBC(VQW.sub(1),p_out,boundaries,4)
bcp5 = DirichletBC(VQW.sub(1),p_out,boundaries,5)

bcp = [bcp6,bcp5]

# SOLID

# MESH DISPLACEMENT

bcU1 = DirichletBC(VQW.sub(2),noslip,boundaries,1)  # walls
bcU2 = DirichletBC(VQW.sub(2),noslip,boundaries,2)  # interface
bcU3 = DirichletBC(VQW.sub(2),noslip,boundaries,3)  # fluid inlet
bcU4 = DirichletBC(VQW.sub(2),noslip,boundaries,4)  # fluid outlet
bcU5 = DirichletBC(VQW.sub(2),noslip,boundaries,5)  # elastic outlet
bcU6 = DirichletBC(VQW.sub(2),noslip,boundaries,6)  # elastic inlet

bcU = [bcU1,bcU5,bcU6]

bcs = [bu for bu in bcu] + [bU for bU in bcU]

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

eta = Function(V)
#u0 = Function(V,'final_u.xml')
u0 = Function(V)
u1 = Function(V)
d1 = Function(V)
UPR = Function(VQW)
n = FacetNormal(mesh)
k = Constant(dt)


nu = 0.658
rho_f = 1.0*1e-3
mu_f = rho_f*nu

# SOLID
E = 5*10**3
Pr = 0.479
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))

rho_p = Constant(1.75e-3)


kappa = Constant(1.4*10**(-15)*(10**9))



phi = 0.2    # porosity
K_perm = kappa/mu_f

rho_s = Constant(rho_p - rho_f*phi/(1-phi))


def sigma_dev(d):
	return 2*mu_s*sym(grad(d)) + lamda*tr(sym(grad(d)))*Identity(2)

def sigma_f(u):
	return 2*mu_f*sym(grad(u))



aW = rho_p/k*inner(d,v)*dx(1,subdomain_data=SD) + \
	rho_f/k*inner(u,v)*dx(1,subdomain_data=SD) - \
	inner(p,div(v))*dx(1,subdomain_data=SD) + \
	k*inner(grad(v),sigma_dev(d))*dx(1,subdomain_data=SD)

LW = rho_f/k*inner(d1,v)*dx(1,subdomain_data=SD) + \
	rho_f/k*inner(u1,v)*dx(1,subdomain_data=SD) - \
	inner(grad(v),sigma_dev(eta))*dx(1,subdomain_data=SD)


aV = rho_f/k*inner(d,w)*dx(1,subdomain_data=SD) + \
	(rho_f/(k*phi)+1./K_perm)*inner(u,w)*dx(1,subdomain_data=SD) - \
	inner(p,div(w))*dx(1,subdomain_data=SD)


LV = rho_f/(k*phi)*inner(u1,w)*dx(1,subdomain_data=SD) + \
	rho_f/k*inner(d1,w)*dx(1,subdomain_data=SD)

aQ = -inner(div(u)+div(d),q)*dx(1,subdomain_data=SD)
LQ = Constant(0)*q*dx(1,subdomain_data=SD)


aS = aV + aW + aQ
LS = LV + LW + LQ

# FLUID

aFW = rho_f/k*inner(u,v)*dx(0,subdomain_data=SD) + \
	rho_f*inner(grad(u0)*(u-d),v)*dx(0,subdomain_data=SD) - \
	 inner(p,div(v))*dx(0,subdomain_data=SD) + \
	mu_f*inner(grad(v),grad(u))*dx(0,subdomain_data=SD)
LFW = rho_f/k*inner(u1,v)*dx(0,subdomain_data=SD)

aFV = inner(grad(d),grad(w))*dx(0,subdomain_data=SD)
LFV = inner(noslip,w)*dx(0,subdomain_data=SD)

aFQ = -inner(div(u),q)*dx(0,subdomain_data=SD)
LFQ = Constant(0)*q*dx(0,subdomain_data=SD)

aF = aFV + aFW + aFQ
LF = LFV + LFW + LFQ

a = aS+aF
L = LS+LF

count = 0
t = dt
while t < T + DOLFIN_EPS:
	p_in.t = t
	p_out.t = t
	
	b = assemble(L)
	eps = 10
	k_iter = 0
	max_iter = 10
	while eps > 1E-7 and k_iter < max_iter:
	    A = assemble(a)
	    A.ident_zeros()
	    [bc.apply(A,b) for bc in bcs]
	    solve(A,UPR.vector(),b,'lu')
	    u_,p_,d_ = UPR.split(True)
	    eps = errornorm(u_,u0,degree_rise=3)
	    k_iter += 1
	    print 'k: ',k_iter, 'error: %.3e' %eps
	    u0.assign(u_)
	
	ufile << u_
	pfile << p_
	dfile << d_
	tfile << eta
	d1.assign(d_)
	d_.vector()[:] *= float(k)
	eta.vector()[:] += d_.vector()[:]
	mesh.move(d_)
	mesh.bounding_box_tree().build(mesh)
        print 't=%.4f'%t
	t += dt

	count += 1





folder = 'initial_data_t_%d'% int(t*10**6)
File('%s/u.xml'%folder) << u_
File('%s/p.xml'%folder) << p_
File('%s/d.xml'%folder) << d_
    
