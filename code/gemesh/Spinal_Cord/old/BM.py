from domains_and_boundaries import *
from pylab import zeros, where, linspace, ones, show, array, loadtxt
import pylab as plt
import scipy.interpolate as Spline
import sys
set_log_active(False)
parameters['allow_extrapolation']=True

'''
RL = sys.argv[1]
EkPa = sys.argv[2]

ref_level = int(RL)
E = Constant(float(EkPa))
'''

ref_level = 1
E = 5.0

print 'ref: ',RL
print 'E: ',EkPa

ufile = File("RESULTS/ref%s_E%s/velocity.pvd"%(RL,EkPa)) # xdmf
pfile = File("RESULTS/ref%s_E%s/pressure.pvd"%(RL,EkPa))
dfile = File("RESULTS/ref%s_E%s/dU.pvd"%(RL,EkPa))
tfile = File("RESULTS/ref%s_E%s/U.pvd"%(RL,EkPa))





Pr = 0.479


rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)

rho_s = Constant(1.75*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))


Nx = 18*ref_level
Ny = 30*ref_level
P1 = Point(-9,0)
P2 = Point(9,60)
mesh = RectangleMesh(P1,P2,Nx,Ny)
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
T = 0.01
# TEST AND TRIALFUNCTIONS

V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = VectorFunctionSpace(mesh,'CG', 1)
VPW = MixedFunctionSpace([V,P,W])
print VPW.dim()
v,p,w = TrialFunctions(VPW)
phi,eta,psi = TestFunctions(VPW)

# PHYSICAL PARAMETERS
h = mesh.hmin()


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0.0,0.0))


C1_Brucker = 10*array([2.26, 1.53, 0.18, 0.19, 0.2, 0.05, -0.39, -0.69, -0.93, -0.66, -0.63, -0.7, -0.23, 2.51])      # 12:00

C1_Brucker2 = 10*array([-0.33, -0.35, -0.51, 0.99, 1.27, 0.83, 0.71, 0.67, -0.15, -0.71, -0.05, -0.21, -0.43, -0.62]) # 6:00


t_Brucker = 1e-3*array([10, 73, 136, 199, 262, 325, 388, 451, 514, 577, 640, 703, 766, 829])
C = C1_Brucker[5:]
C2 = C1_Brucker[:5]
C1_Brucker = plt.append(C,C2)

Cop = C1_Brucker2[5:]
Cop2 = C1_Brucker2[:5]
C1_Brucker2 = plt.append(Cop,Cop2)

Bertram09_Pressure = array([0, -100, 0, 0])
B09time = array([0, 2.5e-3, 5e-3, 0.829])

class MyExpression0(Expression):
	def __init__(self,t_Brucker,C1_Brucker,t,amp=1):
		self.t = t
		self.t_Brucker = t_Brucker
		self.C1_Brucker = C1_Brucker
		self.amp = amp

	def eval(self,values,x):
		t = self.t
		t_Brucker = self.t_Brucker
		C1_Brucker = self.C1_Brucker
		while t > self.t_Brucker[-1]:
			t -= self.t_Brucker[-1]
		tval = t_Brucker
		yval = C1_Brucker
		idx = plt.find(t_Brucker >= t)[0] - 1
		values[0] = self.amp*Spline.UnivariateSpline(tval,yval)(t)

class MyExpression1(Expression):
	def __init__(self,t_Brucker,C1_Brucker,t,amp=1):
		self.t = t
		self.t_Brucker = t_Brucker
		self.C1_Brucker = C1_Brucker
		self.amp = amp

	def eval(self,values,x):
		t = self.t
		t_Brucker = self.t_Brucker
		C1_Brucker = self.C1_Brucker
		while t > self.t_Brucker[-1]:
			t -= self.t_Brucker[-1]
		tval = t_Brucker
		yval = C1_Brucker
		idx = plt.find(t_Brucker >= t)[0] - 1
		values[1] = self.amp*Spline.UnivariateSpline(tval,yval)(t)
		values[0] = 0
	def value_shape(self):
		return (2,)


[T_Eide, P_Eide] = loadtxt('Eide_normalized.txt').transpose()
P_Eide *= -1#33*6./50    # Convert to Pa and make downwards

t_Erika = linspace(0,1.1,23)
P_Erika = array([-0.011,-0.03,-0.02, 0.002,-0.001, -0.002,-0.003, -0.004,0.001, 0.002,0.003, 0.003, 0.004, 0.004,0.003, 0.004,0.006, 0.04,0.045, 0.01, -0.01,-0.01,-0.01])
P_Erika *= 133*6
#pressure = MyExpression0(t_Erika,P_Erika,0.0,amp=1)

pressure = Expression(('amp*sin(2*pi*t)'),t=0,amp=1)


bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Fluid out
bcv5 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Solid out
bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Interface
bcv7 = DirichletBC(VPW.sub(0),noslip,boundaries,7) # Fluid walls


bcv = [bcv3, bcv5, bcv7] # don't use bcv6 for FSI


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

#v1 = Expression(('0','6*(x[0]+L)*(x[0]+9)'),L=5)
#v1 = interpolate(v1,V)
v1 = Function(V)#,'initial_data/u.xml')
U = Function(W)
U0 = Function(W)


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


epsilon = 1e8


U_ = U + k*v      # discretize in time here

D = inner(sigma_dev(U_),grad(phi))*dx_s

aMS = rho_s/k*inner(v,phi)*dx_s \
	+ rho_s*inner(grad(v0)*v,phi)*dx_s \
	+ lhs(D)# inner(sigma_dev(k*v),grad(phi))*dx_s

LMS = rho_s/k*inner(v1,phi)*dx_s \
	+ rhs(D)#inner(sigma_dev(U),grad(phi))*dx_s
'''
aMS = rho_s/k*inner(v,phi)*dx_s \
	+ rho_s*inner(grad(v0)*v,phi)*dx_s \
	+ k*inner(sigma_dev(v),grad(phi))*dx_s

LMS = rho_s/k*inner(v1,phi)*dx_s - \
	inner(sigma_dev(U),grad(phi))*dx_s# - \
	#inner(g,phi)*dx_s
'''
aDS = epsilon*inner(v,psi)*dx_s - epsilon*inner(w,psi)*dx_s



aS = aMS + aDS
LS = LMS

# FLUID
penalty = 0.05*mesh.hmin()


aMF = rho_f/k*inner(v,phi)*dx_f \
	+ rho_f*inner(grad(v0)*(v-w),phi)*dx_f \
	+ inner(p,div(phi))*dx_f + \
	2*mu_f*inner(sym(grad(v)),grad(phi))*dx_f \
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

aCF = -inner(div(v),eta)*dx_f

aF = aMF + aDF + aCF
LF = LMF + LDF

a = aS+aF
L = LS+LF

'''
if initial_files:
	t = t_folder*1e-6
	with open('%s/numbers.txt'%folder,'r') as infile:
		xA,yA,coord = infile.readline().split()
		xA,yA,coord = float(xA),float(yA),int(coord)	

else:
'''
t = dt
count = 0



while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):

	if t < 2.0:
		pressure.amp = 5*t
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
	if count%10==0:
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
	U0.assign(U)
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
