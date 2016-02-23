from domains_and_boundaries import *
from Brucker_Inlet_Vel import *
import sys
#set_log_active(False)
parameters['allow_extrapolation']=True

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
	mesh = Mesh('meshes/wide_syrinx.xml')#cord_w_csc_mm.xml')
	#mesh = refine(mesh)
	#mesh = refine(mesh)
	p1 = Point(-9,0)
	p2 = Point(-5,60)
	#mesh = RectangleMesh(p1,p2,8,60)
	LEN = len(mesh.coordinates())
	print 'len(mesh.coord)', LEN
	SD = CellFunction('size_t', mesh)
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


dt = 0.001   # use 0.0003 for oscillations
T = 3.0
# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = FunctionSpace(mesh,'CG', 2)
VPW = MixedFunctionSpace([V,P,V])

VPW_ = Function(VPW)

v,p,w = split(VPW_)
phi,eta,psi = TestFunctions(VPW)

# PHYSICAL PARAMETERS
h = mesh.hmin()


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0.0,0.0))


vf_l = MyExpression0(t_Brucker,C1_Brucker,0.0)
vf_r = MyExpression0(t_Brucker,C1_Brucker,0.0)

slip = Constant((0.1,0))

bcv1 = DirichletBC(VPW.sub(0),vf_l,boundaries,1)     # Fluid in_l
bcv2 = DirichletBC(VPW.sub(0),vf_r,boundaries,2)	# Fluid in_r
bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Fluid out
bcv5 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Solid out
bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Interface
bcv7 = DirichletBC(VPW.sub(0),noslip,boundaries,7) # Fluid walls


bcv = [bcv1, bcv2, bcv3, bcv5, bcv7] # don't use bcv6 for FSI

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
	U = Function(V,'%s/U.xml'%folder)
else:
	VPW1 = Function(VPW)
	v1,p1,w1 = split(VPW1)
	U = Function(V)

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

#AV = epsilon*inner(grad(u),grad(w))*dx_s - epsilon*inner(grad(d),grad(w))*dx_s
aDS = epsilon*inner(v,psi)*dx_s - epsilon*inner(w,psi)*dx_s


aCS = epsilon*p*eta*dx_s

aS = aMS + aDS + aCS
LS = LMS

# FLUID
aMF = rho_f/k*inner(v,phi)*dx_f + \
	rho_f*inner(grad(v)*(v-w),phi)*dx_f - \
	 inner(p,div(phi))*dx_f + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx_f

LMF = rho_f/k*inner(v1,phi)*dx_f

aDF = k*inner(grad(w),grad(psi))*dx_f
LDF = -inner(grad(U),grad(psi))*dx_f

aCF = -inner(div(v),eta)*dx_f
LFQ = Constant(0)*eta*dx_f

aF = aMF + aCF + aDF
LF = LMF + LDF + LFQ

a = aS+aF
L = LS+LF

F = a+L

if initial_files:
	t = t_folder*1e-6
	with open('%s/numbers.txt'%folder,'r') as infile:
		xA,yA,coord = infile.readline().split()
		xA,yA,coord = float(xA),float(yA),int(coord)	

else:
	t = dt
	xA = where(mesh.coordinates()[:,0] == 0.6)
	yA = where(mesh.coordinates()[:,1] == 0.2)

	n1 = as_vector((1.0,0))
	n2 = as_vector((0,1.0))
	nx = dot(n1,n)
	ny = dot(n2,n)
	nt = as_vector((ny,-nx))

count = 0


change_vf = True
Fd = plt.zeros(1e6)
Fl = plt.zeros(1e6)
Fd[-1], Fl[-1] = 1,1
FdC, FlC = 1,1
ydisp = []
xdisp = []
time = []
n1 = as_vector((1.0,0))
n2 = as_vector((0,1.0))
nx = dot(n1,n)
ny = dot(n2,n)
nt = as_vector((ny,-nx))


while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):
	vf_l.t=t
	vf_r.t=t
	solve(F==0,VPW_,bcs=bcv+bcw)
	v1,p1,w1 = split(VPW_)#.split(True)
	v_,p_,w_ = VPW_.split(True)
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
	VPW1.assign(VPW_)
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
