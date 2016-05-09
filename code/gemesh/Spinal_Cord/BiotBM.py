from domains_and_boundaries import *
from pylab import zeros, where, linspace, ones, show, array
import pylab as plt
import sys
set_log_active(False)
parameters['allow_extrapolation']=True

ufile = File("results_mono/velocity.pvd") # xdmf
pfile = File("results_mono/pressure.pvd")
dfile = File("results_mono/dU.pvd")
tfile = File("results_mono/U.pvd")
qfile = File("results_mono/q.pvd")

Pr = 0.479

RL = '3'
EkPa = '5000'
ref_level = int(RL)
E = Constant(float(EkPa))

rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)

rho_s = Constant(1.75*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))


initial_files = False
if initial_files:
	t_folder = 3000098
	folder = 'initial_data_t_%d' % t_folder
	mesh = Mesh('%s/mesh.xml'%folder)
	SD = MeshFunction('size_t',mesh,'%s/mesh_func.xml'%folder)
	boundaries = MeshFunction('size_t',mesh,'%s/facet_func.xml'%folder)

else:
	P0 = Point(-9,0)
	P1 = Point(9,60)
	mesh = RectangleMesh(P0,P1,ref_level*18,ref_level*30)
	#mesh = refine(mesh)
	#mesh = refine(mesh)
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
	CSC_bnd().mark(boundaries,6)
# Pr = 0.479
# x max for U is 3.17*10^-5 for E = 5*10^6 and rho = 0.001*rho_f at time 0.4 (timestep 40 for dt = 0.01)
# x max for U is 3.17*10^-5 for E = 5*10^6 and rho = 2*rho_f at time 0.4 (timestep 40 for dt = 0.01)
# x max for U is 2.99*10^-5 for E = 5*10^6 and rho = 1000*rho_f at time 0.4 (timestep 40 for dt = 0.01)


# E = 5*10^6, rho = 2*rho_f, time = 0.4, dt = 0.01, (timestep 40)
# x max for U is 3.09*10^-5 for Pr = 0.4999
# x max for U is 3.45*10^-5 for Pr = 0.4

#E = 16*10^3, rho = 2*rho_f, time = 0.4, dt = 0.01, (timestep 40)
# x max for U is 4.69^-3 for Pr = 0.4999





dt = 0.002   # use 0.0003 for oscillations
T = 10
# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = VectorFunctionSpace(mesh,'CG', 2)
VPW = MixedFunctionSpace([V,P,W])

v,p,w = TrialFunctions(VPW)
phi,eta,psi = TestFunctions(VPW)

# PHYSICAL PARAMETERS
h = mesh.hmin()


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0.0,0.0))




vf_l = Constant((0,-0.1))
vf_r = Constant((0,-0.1))
pressure = Expression(('amp*sin(2*pi*t)'),t=0,amp=1)

slip = Constant((0.1,0))

bcv1 = DirichletBC(VPW.sub(0),vf_l,boundaries,1)     # Fluid in_l
bcv2 = DirichletBC(VPW.sub(0),vf_r,boundaries,2)	# Fluid in_r
bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Fluid out
bcv5 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Solid out
bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Interface
bcv7 = DirichletBC(VPW.sub(0),noslip,boundaries,7) # Fluid walls

bcp = DirichletBC(VPW.sub(1),Constant(1),boundaries,1)
bcp2 = DirichletBC(VPW.sub(1),Constant(1),boundaries,2)


bcv = [bcv3,bcv5,bcv7] # don't use bcv6 for FSI

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
v0 = Function(V)
w0 = Function(V)

if initial_files:
	v1 = Function(V,'%s/u.xml'%folder)
	U = Function(W,'%s/U.xml'%folder)
	
	#plot(u1)
else:
	v1 = Function(V)#,'initial_data/u.xml')
	U = Function(W)
	w1 = Function(W)
	q = Function(V)


VPW_ = Function(VPW)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))
n = FacetNormal(mesh)
g = Constant((0,2.0))


dS = Measure('dS')[boundaries]
dx = Measure('dx')[SD]
ds = Measure('ds')[boundaries]

dx_f = dx(0,subdomain_data=SD)
dx_s = dx(1,subdomain_data=SD)

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
penalty = 0.05*mesh.hmin()
gamma = Constant(1.0)
poro = Constant(0.2)
K_D = K_perm/mu_f
rho_p = Constant(rho_s*(1-poro) + rho_f*poro)


aMS = rho_p/k*inner(w,phi)*dx_s \
	+ rho_p/k*inner(grad(w0)*w,phi)*dx_s \
	+ k*inner(sigma_dev(w),grad(phi))*dx_s \
	- inner(p,div(phi))*dx_s

LMS = rho_p/k*inner(w1,phi)*dx_s - \
	inner(sigma_dev(U),grad(phi))*dx_s
	# - \
	#inner(pressure*n,phi)*ds(3) + \
	#inner(pressure*n,phi)*ds(5)


#total velocity flux not experienced velocity

# filtration vel: q = -K*grad(p) - K*rho_f*dw/dt 
# v = w_s + q
# v = w_s - K*grad(p) - K*rho_f*dw/dt
# div (v) = 0

aDS = epsilon*inner(v,psi)*dx_s \
	- epsilon*inner(w,psi)*dx_s \
	+ epsilon*K_D*inner(grad(p),psi)*dx_s \
	+ epsilon*rho_f*K_D/k*inner(w,psi)*dx_s \
	+ epsilon*rho_f*K_D*inner(grad(w0)*w,psi)*dx_s

LDS = epsilon*rho_f*K_D/k*inner(w1,psi)*dx_s


aCS = -inner(div(v),eta)*dx_s

aS = aMS + aDS + aCS
LS = LMS + LDS

# FLUID

aMF = rho_f/k*inner(v,phi)*dx_f + \
	rho_f*inner(grad(v0)*(v-w),phi)*dx_f - \
	 inner(p,div(phi))*dx_f + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx_f \
	- mu_f*inner(grad(v).T*n,phi)*ds(1) \
	- mu_f*inner(grad(v).T*n,phi)*ds(2) \
	- mu_f*inner(grad(v).T*n,phi)*ds(4) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(1) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(2) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(4) #\
	#- gamma/sqrt(K_D)*inner((v-w)-dot(v-w,n)*n,(phi-psi)-dot(phi-psi,n)*n)*ds(6)

LMF = rho_f/k*inner(v1,phi)*dx_f - \
	  inner(pressure*n,phi)*ds(1) - \
	  inner(pressure*n,phi)*ds(2) + \
	  inner(pressure*n,phi)*ds(4)

aDF = k*inner(grad(w),grad(psi))*dx_f \
	+ k*inner(grad(w('-'))*n('-'),psi('-'))*dS(6)
LDF = -inner(grad(U),grad(psi))*dx_f \
	+ inner(grad(U('-'))*n('-'),psi('-'))*dS(6)
aCF = -inner(div(v),eta)*dx_f

aF = aMF + aDF + aCF
LF = LMF + LDF

a = aS+aF
L = LS+LF



count = 0


#accumulated = 0
#flowfile = open('out.txt','w')

t=dt

while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):
	#vf_l.t=t
	#vf_r.t=t	
	if t < 1.0:
		pressure.amp = 10*t
	pressure.t=t
	b = assemble(L)
	eps = 10
	k_iter = 0
	max_iter = 5
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
	    w0.assign(w_)
	q.vector()[:] = v_.vector()[:]
	q.vector()[:] -= w_.vector()[:]
	if count%10==0:
		qfile << q
		ufile << v_
		pfile << p_
		dfile << w_
		tfile << U
	w_.vector()[:] *= float(k)
	U.vector()[:] += w_.vector()[:]
	mesh.move(w_)
	mesh.bounding_box_tree().build(mesh)
	#print '%10.4e %5.4e' %(XD, YD)
	# Move to next time step
	v1.assign(v_)
	w1.assign(w_)
	
	#flow = assemble(dot(q('+'),n('+'))*dS(6))
	#accumulated += flow
	
        print 't=%.4f'%t#   flow=%g'%(t,flow)
	t += dt
	#print '%10d %10.4e %5.4e' %(LEN,xA-mesh.coordinates()[coord,0], yA-mesh.coordinates()[coord,1])
	count += 1
	#print 'accumulated = %g'%accumulated
	#flowfile.write('%g   %g'%(flow,accumulated))

#flowfile.close()

#measure max q by extraction region and fit to temporal range?


'''
folder = 'initial_data_t_%d'% int(t*10**6)
File('%s/u.xml'%folder) << v_
File('%s/p.xml'%folder) << p_
File('%s/d.xml'%folder) << w_
File('%s/U.xml'%folder) << U
File('%s/mesh.xml'%folder) << mesh
File('%s/mesh_func.xml'%folder) << SD
File('%s/facet_func.xml'%folder) << boundaries   
with open('%s/numbers.txt'%folder,'w') as outfile:
	outfile.write('%g %g %g' %(xA,yA,coord))



from pylab import plot,show, figure
plot(time,xdisp)
figure()
plot(time,ydisp)
figure()
plot(time,Fd[:count])
figure()
plot(time,Fl[:count])
show()
'''
