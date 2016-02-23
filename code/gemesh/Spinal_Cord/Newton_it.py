from domains_and_boundaries import *
from pylab import zeros, where
import pylab as plt
import sys

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
	mesh = Mesh('meshes/cord_w_csc_mm.xml')
	#mesh = refine(mesh)
	#mesh = refine(mesh)
	LEN = len(mesh.coordinates())
	print 'len(mesh.coord)', LEN
	SD = MeshFunction('size_t', mesh, mesh.topology().dim())
	SD.set_all(0)
	#Solid().mark(SD,1)
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


dt = 0.0003   # use 0.0003 for oscillations
T = 10*dt
# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = FunctionSpace(mesh,'CG', 2)
VPW = MixedFunctionSpace([V,P,V])

vpw = TrialFunction(VPW)
pep = TestFunction(VPW)

v,p,w = split(vpw)
phi,eta,psi = split(pep)

# PHYSICAL PARAMETERS
h = mesh.hmin()


# INITIAL AND BOUNDARY CONDITIONS

# FLUID
noslip = Constant((0.0,0.0))
vf_l = Expression(('0.0','10*std::exp(-pow(t-0.5,2)/0.1)*(x[0]+x0)*(x1+x[0])'),x0=x0,x1=x1,t=dt)
vf_r = Expression(('0.0','-10*std::exp(-pow(t-0.5,2)/0.1)*(x[0]-x0)*(x1-x[0])'),x0=x0,x1=x1,t=dt)

slip = Constant((0.1,0))

bcv1 = DirichletBC(VPW.sub(0),vf_l,boundaries,1)     # Fluid in_l
bcv2 = DirichletBC(VPW.sub(0),vf_r,boundaries,2)	# Fluid in_r
bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid in
bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Fluid out
bcv5 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Solid out
bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Interface
bcv7 = DirichletBC(VPW.sub(0),noslip,boundaries,7) # Fluid walls


bcv = [bcv1, bcv2, bcv3, bcv7]

# SOLID

# MESH DISPLACEMENT

bcw1 = DirichletBC(VPW.sub(2),noslip,boundaries,1)  # Fluid in_l
bcw2 = DirichletBC(VPW.sub(2),noslip,boundaries,2)  # Fluid in_r
bcw3 = DirichletBC(VPW.sub(2),noslip,boundaries,3)  # Solid in
bcw4 = DirichletBC(VPW.sub(2),noslip,boundaries,3)  # Fluid out
bcw5 = DirichletBC(VPW.sub(2),noslip,boundaries,5)  # Solid out
bcw6 = DirichletBC(VPW.sub(2),noslip,boundaries,6)  # Interface
bcw7 = DirichletBC(VPW.sub(2),noslip,boundaries,7) # Fluid walls

bcw = [bcw1,bcw2,bcw3,bcw5,bcw7]

# CREATE FUNCTIONS
v0 = Function(V)
if initial_files:
	v1 = Function(V,'%s/u.xml'%folder)
	U = Function(V,'%s/U.xml'%folder)
	#plot(u1)
else:
	v1 = Function(V)#,'initial_data/u.xml')
	U = Function(V)


VPW_ = Function(VPW)
dVPW_ = Function(VPW)

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
'''
aMS = rho_f/k*inner(v,phi)*dx(1,subdomain_data=SD) + \
	rho_f*inner(grad(v)*(v-w),phi)*dx(1,subdomain_data=SD) - \
	 inner(p,div(phi))*dx(1,subdomain_data=SD) + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx(1,subdomain_data=SD)
'''
LCS = inner(div(v),eta)*dx(1,subdomain_data=SD)

aS = aMS + aDS
LS = LMS + LCS

# FLUID
'''
aFW = rho_f/k*inner(u,v)*dx(0,subdomain_data=SD) + \
	rho_f*inner(grad(u0)*(u-d),v)*dx(0,subdomain_data=SD) - \
	 inner(p,div(v))*dx(0,subdomain_data=SD) + \
	mu_f*inner(grad(v),grad(u))*dx(0,subdomain_data=SD)
'''
aMF = rho_f/k*inner(v,phi)*dx(0,subdomain_data=SD) + \
	rho_f*inner(grad(v)*(v-w),phi)*dx(0,subdomain_data=SD) - \
	 inner(p,div(phi))*dx(0,subdomain_data=SD) + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx(0,subdomain_data=SD)

LMF = rho_f/k*inner(v1,phi)*dx(0,subdomain_data=SD)

aDF = k*inner(grad(w),grad(psi))*dx(0,subdomain_data=SD)
LDF = -inner(grad(U),grad(psi))*dx(0,subdomain_data=SD)

aCF = -inner(div(v),eta)*dx(0,subdomain_data=SD)
#LFQ = Constant(0)*eta*dx(0,subdomain_data=SD)

aF = aMF + aCF + aDF
LF = LMF + LDF 

a = aS+aF
L = LS+LF

F = a-L
F = action(F,VPW_)
J = derivative(F,VPW_,vpw)

problem = NonlinearVariationalProblem(F,VPW_, bcv+bcw,J)
class MySolver(NonlinearVariationalSolver):
	def __init__(self,problem):
		NonlinearVariationalSolver.__init__(self,problem)

	def solve(self):
		print self.__init__.func_code.co_varnames
		NonlinearVariationalSolver.solve(self)


solver = MySolver(problem)
prm = solver.parameters
#prm['newton_solver']['absolute_tolerance'] = 1E1
#prm['newton_solver']['relative_tolerance'] = 1E1
solver.solve()
v_,p_,w_ = VPW_.split(True)
plot(p_)
interactive()
print 'OK'
sys.exit()


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
	xA = plt.where(abs(mesh.coordinates()[:,0] - 0.6) < h)
	yA = plt.where(abs(mesh.coordinates()[:,1] - 0.2) < h)
	for idx in xA[0]:
		if idx in yA[0]:
		    coord = idx
	xA = float(mesh.coordinates()[coord,0])
	yA = float(mesh.coordinates()[coord,1])
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
	if t <= 2.0:
		vf_l.t=t
		vf_r.t=t
	else:
		vf.t = 2.0
	b = assemble(L)
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
	'''
	D = interpolate(w_,V)
	uf = interpolate(v_,V)
	pf = interpolate(p_,P)
	ut = dot(uf,nt)
	DragSolid = -assemble((rho_f*nu*dot(n('-'),grad(ut('-')))*ny('-') - pf('-')*nx('-'))*dS(5))
	LiftSolid = assemble((rho_f*nu*dot(n('-'),grad(ut('-')))*nx('-') + pf('-')*ny('-'))*dS(5))
	ut = dot(uf,nt)
	Drag = -assemble((rho_f*nu*dot(n,grad(ut))*ny - pf*nx)*ds(6))
	Lift = assemble((rho_f*nu*dot(n,grad(ut))*nx + pf*ny)*ds(6))
	print 'Drag: ', Drag + DragSolid
	print 'Lift: ', Lift + LiftSolid
	#bcu[-1] = DirichletBC(VQW.sub(0),d_,boundaries,5) # interface
	Fd[count] = DragSolid + Drag
	Fl[count] = LiftSolid + Lift
	FdC = 100*(Fd[count]/Fd[count-1]-1)
	FlC = 100*(Fl[count]/Fl[count-1]-1)
	#print 'percent change: Fd: %.3f       Fl: %.3f' %(FdC, FlC)
	YD = mesh.coordinates()[coord,1]-yA
	XD = mesh.coordinates()[coord,0]-xA
	ydisp.append(YD)
	xdisp.append(XD)
	time.append(t)
	'''
	if count%1==0:
		ufile << v_
		pfile << p_
		dfile << w_
		tfile << U
	#w_.vector()[:] *= float(k)
	U.vector()[:] += float(k)*w_.vector()[:]
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
