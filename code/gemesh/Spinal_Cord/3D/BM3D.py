from domains_3D import *
from pylab import zeros, where, linspace, ones, show, array, loadtxt
import pylab as plt
import scipy.interpolate as Spline
import sys
set_log_active(False)
parameters['allow_extrapolation']=True


RL = '1'
EkPa = '16000'
ref_level = int(RL)
E = Constant(float(EkPa))

print 'ref: ',RL
print 'E: ',EkPa






Pr = 0.479


rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)

rho_s = Constant(1.75*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))



mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "drosdal3D.h5", "r")
hdf.read(mesh, "/mesh", False)
SD = CellFunction("size_t", mesh)
hdf.read(SD, "/subdomains")
boundaries = FacetFunction("size_t", mesh)
hdf.read(boundaries, "/boundaries")

'''
mesh = Mesh('drosdal3D.xml')
LEN = len(mesh.coordinates())
print 'len(mesh.coord)', LEN

SD = MeshFunction('size_t', mesh, 'drosdal3D_physical_region.xml')
# DEFINING BOUNDARIES
boundaries = MeshFunction("size_t",mesh,'drosdal3D_facet_region.xml')
'''

dt = 0.002   # use 0.0003 for oscillations
T = 10
# TEST AND TRIALFUNCTIONS

V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = VectorFunctionSpace(mesh,'CG', 1)
VPW = MixedFunctionSpace([V,P,W])
print VPW.dim()
v,p,w = TrialFunctions(VPW)
phi,eta,psi = TestFunctions(VPW)


# PHYSICAL PARAMETERS
hmesh = mesh.hmin()


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0.0,0.0,0.0))

pressure = Expression(('amp*sin(2*pi*t)'),t=0,amp=1)

# fluid out
# fluid in
# solid out 
# solid in
# fluid walls
# interface

bcv2 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid out
bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Fluid walls

bcv = [bcv2, bcv3, bcv4] # don't use bcv6 for FSI


# SOLID

# MESH DISPLACEMENT
bcw1 = DirichletBC(VPW.sub(2),noslip,boundaries,1) # Fluid out
bcw2 = DirichletBC(VPW.sub(2),noslip,boundaries,2) # Fluid in
bcw3 = DirichletBC(VPW.sub(2),noslip,boundaries,3) # Solid out
bcw4 = DirichletBC(VPW.sub(2),noslip,boundaries,4) # Solid in
bcw5 = DirichletBC(VPW.sub(2),noslip,boundaries,5) # Fluid walls


bcw = [bcw1,bcw2,bcw3,bcw4,bcw5]

# CREATE FUNCTIONS
v0 = Function(V)

#v1 = Expression(('0','6*(x[0]+L)*(x[0]+9)'),L=5)
#v1 = interpolate(v1,V)
v1 = Function(V)#,'initial_data/u.xml')
U = Function(W)


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


def sigma_dev(U):
	return 2*mu_s*sym(grad(U)) + lamda*tr(sym(grad(U)))*Identity(3)


epsilon = 1e8


U_ = U + k*v      # discretize in time here

D = inner(sigma_dev(U_),grad(phi))*dx_s

aMS = rho_s/k*inner(v,phi)*dx_s \
	+ rho_s*inner(grad(v0)*v,phi)*dx_s \
	+ lhs(D)# inner(sigma_dev(k*v),grad(phi))*dx_s

LMS = rho_s/k*inner(v1,phi)*dx_s \
	+ rhs(D)#inner(sigma_dev(U),grad(phi))*dx_s

aDS = epsilon*inner(v,psi)*dx_s - epsilon*inner(w,psi)*dx_s



aS = aMS + aDS
LS = LMS

# FLUID
penalty = 0.01*hmesh


aMF = rho_f/k*inner(v,phi)*dx_f \
	+ rho_f*inner(grad(v0)*(v-w),phi)*dx_f \
	- inner(p,div(phi))*dx_f + \
	2*mu_f*inner(sym(grad(v)),grad(phi))*dx_f \
	- mu_f*inner(grad(v).T*n,phi)*ds(1) \
	- mu_f*inner(grad(v).T*n,phi)*ds(0) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(1) \
	+ penalty**-2*(inner(v,phi)-inner(dot(v,n),dot(phi,n)))*ds(4)

LMF = rho_f/k*inner(v1,phi)*dx_f - \
	inner(pressure*n,phi)*ds(1) + \
	inner(pressure*n,phi)*ds(0)

aDF = k*inner(grad(w),grad(psi))*dx_f
LDF = -inner(grad(U),grad(psi))*dx_f

aCF = -inner(div(v),eta)*dx_f

aF = aMF + aDF + aCF
LF = LMF + LDF

a = aS+aF
L = LS+LF

t = dt
count = 0

ufile = File("RESULTS/ref%s_E%s/velocity.pvd"%(RL,EkPa)) # xdmf
pfile = File("RESULTS/ref%s_E%s/pressure.pvd"%(RL,EkPa))
dfile = File("RESULTS/ref%s_E%s/dU.pvd"%(RL,EkPa))
tfile = File("RESULTS/ref%s_E%s/U.pvd"%(RL,EkPa))


while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):

	if t < 1.0:
		pressure.amp = 10*t
	pressure.t = t
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
	if count%5==0:
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
        print 't=%.4f'%t
	t += dt
	#print '%10d %10.4e %5.4e' %(LEN,xA-mesh.coordinates()[coord,0], yA-mesh.coordinates()[coord,1])
	count += 1
	
	
	#compute_norm_at_vertices(v_,norm)
	#print norm.vector().array()




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
