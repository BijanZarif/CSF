from domains_and_boundaries import *
from pylab import zeros, where, linspace, ones, show, array, loadtxt
import pylab as plt
import scipy.interpolate as Spline
import sys
set_log_active(False)
parameters['allow_extrapolation']=True




Pr = 0.479
E = 16

rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)

rho_s = Constant(1.75*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))


ufile = File("results_mono/velocity.pvd") # xdmf
pfile = File("results_mono/pressure.pvd")
dfile = File("results_mono/dU.pvd")
tfile = File("results_mono/U.pvd")

initial_files = False
if initial_files:
	t_folder = 3000098
	folder = 'initial_data_t_%d' % t_folder
	mesh = Mesh('%s/mesh.xml'%folder)
	SD = MeshFunction('size_t',mesh,'%s/mesh_func.xml'%folder)
	boundaries = MeshFunction('size_t',mesh,'%s/facet_func.xml'%folder)

else:
	#mesh = Mesh('meshes/wide_syrinx.xml')#cord_w_csc_mm.xml')
	P1 = Point(-9,0)
	P2 = Point(9,60)
	mesh = RectangleMesh(P1,P2,8,80)
	
	#mesh = refine(mesh)
	#mesh = refine(mesh)
	LEN = len(mesh.coordinates())
	print 'len(mesh.coord)', LEN
	SD = MeshFunction('size_t', mesh, mesh.topology().dim())
	SD.set_all(0)
	Solid().mark(SD,1)
	#CSC().mark(SD,1)

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

# Pr = 0.479
# x max for U is 3.17*10^-5 for E = 5*10^6 and rho = 0.001*rho_f at time 0.4 (timestep 40 for dt = 0.01)
# x max for U is 3.17*10^-5 for E = 5*10^6 and rho = 2*rho_f at time 0.4 (timestep 40 for dt = 0.01)
# x max for U is 2.99*10^-5 for E = 5*10^6 and rho = 1000*rho_f at time 0.4 (timestep 40 for dt = 0.01)


# E = 5*10^6, rho = 2*rho_f, time = 0.4, dt = 0.01, (timestep 40)
# x max for U is 3.09*10^-5 for Pr = 0.4999
# x max for U is 3.45*10^-5 for Pr = 0.4

#E = 16*10^3, rho = 2*rho_f, time = 0.4, dt = 0.01, (timestep 40)
# x max for U is 4.69^-3 for Pr = 0.4999





dt = 0.001   # use 0.0003 for oscillations
T = 10.0
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
vf_l = Expression(('0.0','20*(std::exp(-pow(t-0.8,2)/0.05)+std::exp(-pow(t-1.7,2)/0.05)+std::exp(-pow(t-2.6,2)/0.05))*(x[0]+x0)*(x1+x[0])/4.0'),x0=x0,x1=x1,t=0)
vf_r = Expression(('0.0','-20*(std::exp(-pow(t-0.8,2)/0.05)+std::exp(-pow(t-1.7,2)/0.05)+std::exp(-pow(t-2.6,2)/0.05))*(x[0]-x0)*(x1-x[0])/4.0'),x0=x0,x1=x1,t=0)

vf_l = Expression(('0.0','-80*(std::exp(-pow(t-0.8,2)/0.05)+std::exp(-pow(t-1.7,2)/0.05)+std::exp(-pow(t-2.6,2)/0.05))'),x0=x0,x1=x1,t=0)
vf_r = Expression(('0.0','-80*(std::exp(-pow(t-0.8,2)/0.05)+std::exp(-pow(t-1.7,2)/0.05)+std::exp(-pow(t-2.6,2)/0.05))'),x0=x0,x1=x1,t=0)

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
		#values[1] = -(yval[idx] + (yval[idx+1]-yval[idx])*(t-tval[idx])/(tval[idx+1]-tval[idx]))
		#values[0] = 0
		#values[1] = amp*Spline.UnivariateSpline(tval,yval)(t)
		values[0] = self.amp*Spline.UnivariateSpline(tval,yval)(t)
	#def value_shape(self):
	#	return (2,)

vf_l = MyExpression0(t_Brucker,C1_Brucker,0.0)
vf_r = MyExpression0(t_Brucker,C1_Brucker,0.0)



vf_l = Constant((0,-0.1))
vf_r = Constant((0,-0.1))
[T_Eide, P_Eide] = loadtxt('Eide_normalized.txt').transpose()
P_Eide *= -1#33*6./50    # Convert to Pa and make downwards

t_Erika = linspace(0,1.1,23)
P_Erika = array([-0.011,-0.03,-0.02, 0.002,-0.001, -0.002,-0.003, -0.004,0.001, 0.002,0.003, 0.003, 0.004, 0.004,0.003, 0.004,0.006, 0.04,0.045, 0.01, -0.01,-0.01,-0.01])
P_Erika *= 133*6
#pressure = MyExpression0(t_Erika,P_Erika,0.0,amp=1)

pressure = Expression(('amp*sin(2*pi*t)'),t=0,amp=1)
pressure2 = Expression(('-amp*sin(2*pi*t)'),t=0,amp=1)

slip = Constant((0.1,0))

bcv1 = DirichletBC(VPW.sub(0),vf_l,boundaries,1)     # Fluid in_l
bcv2 = DirichletBC(VPW.sub(0),vf_r,boundaries,2)	# Fluid in_r
bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Fluid out
bcv5 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Solid out
bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Interface
bcv7 = DirichletBC(VPW.sub(0),noslip,boundaries,7) # Fluid walls


bcv = [bcv3, bcv5, bcv7] # don't use bcv6 for FSI

bcp1 = DirichletBC(VPW.sub(1),pressure,boundaries,1)
bcp2 = DirichletBC(VPW.sub(1),pressure,boundaries,2)
bcp4 = DirichletBC(VPW.sub(1),pressure2,boundaries,4)
bcp = [bcp1,bcp2,bcp4]

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
if initial_files:
	v1 = Function(V,'%s/u.xml'%folder)
	U = Function(W,'%s/U.xml'%folder)
	#plot(u1)
else:
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

aMS = rho_s/k*inner(v,phi)*dx_s + \
	k*inner(sigma_dev(v),grad(phi))*dx_s

LMS = rho_s/k*inner(v1,phi)*dx_s - \
	inner(sigma_dev(U),grad(phi))*dx_s# - \
	#inner(g,phi)*dx_s

aDS = epsilon*inner(v,psi)*dx_s - epsilon*inner(w,psi)*dx_s


aCS = -inner(div(v),eta)*dx_s


aS = aMS + aDS# + aCS
LS = LMS

# FLUID

aMF = rho_f/k*inner(v,phi)*dx_f + \
	rho_f*inner(grad(v0)*(v-w),phi)*dx_f - \
	 inner(p,div(phi))*dx_f + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx_f - \
	mu_f*inner(grad(v).T*n,phi)*ds(1) - \
	mu_f*inner(grad(v).T*n,phi)*ds(4) - \
	mu_f*inner(grad(v).T*n,phi)*ds(2)

LMF =rho_f/k*inner(v1,phi)*dx_f - \
	inner(pressure*n,phi)*ds(1) - \
	inner(pressure*n,phi)*ds(2) + \
	inner(pressure*n,phi)*ds(4)

aCF = -inner(div(v),eta)*dx_f
aDF = k*inner(grad(w),grad(psi))*dx_f
LDF = -inner(grad(U),grad(psi))*dx_f

aF = aMF + aCF + aDF
LF = LMF + LDF

a = aF + aS
L = LF + LS
F = a+L

if initial_files:
	t = t_folder*1e-6
	with open('%s/numbers.txt'%folder,'r') as infile:
		xA,yA,coord = infile.readline().split()
		xA,yA,coord = float(xA),float(yA),int(coord)	

else:
	t = dt
	'''
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
	'''
count = 0


change_vf = True
Fd = plt.zeros(1e6)
Fl = plt.zeros(1e6)
Fd[-1], Fl[-1] = 1,1
FdC, FlC = 1,1
ydisp = []
xdisp = []
time = []
'''
n1 = as_vector((1.0,0))
n2 = as_vector((0,1.0))
nx = dot(n1,n)
ny = dot(n2,n)
nt = as_vector((ny,-nx))
'''
while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):
	#vf_l.t=t
	#vf_r.t=t
	#if t < 2.0:
	#	pressure.amp=t/2.0
	pressure.t = 0.25
	b = assemble(L)
	eps = 10
	k_iter = 0
	max_iter = 5
	while eps > 1E-10 and k_iter < max_iter:
	    A = assemble(a)
	    A.ident_zeros()
	    [bc.apply(A,b) for bc in bcv]
	    [bc.apply(A,b) for bc in bcw]
	    #[bc.apply(A,b) for bc in bcp]
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


'''
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
