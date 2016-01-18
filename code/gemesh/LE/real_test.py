from dolfin import *
set_log_active(False)
import sys

h = 6
L = 60

p0 = Point(0,0)
p1 = Point(L,2*h)

ufile = File("results_mono/velocity.pvd")
pfile = File("results_mono/pressure.pvd")
dfile = File("results_mono/dU.pvd")
tfile = File("results_mono/U.pvd")

dt = 0.001
T = 1000*dt

mesh = RectangleMesh(p0,p1,40,8)

eps = 1e-10

class InletS(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[0],0.0)  and x[1] > h - eps
class OutletS(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[0],L) and x[1] > h - eps
class InletF(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[0],0.0)  and x[1] < h + eps
class OutletF(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[0],L) and x[1] < h + eps
class Top(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[1],2*h)
class Bottom(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd and near(x[1],0.0)

class Interface(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],h)

class Fluid(SubDomain):
	def inside(self,x,on_bnd):
		return x[1] < h + eps

# DOMAIN
SD = MeshFunction('uint',mesh,mesh.topology().dim())
SD.set_all(1)
Fluid().mark(SD,0)

bnd = FacetFunction("size_t",mesh)
InletS().mark(bnd,1)
OutletS().mark(bnd,2)
InletF().mark(bnd,3)
OutletF().mark(bnd,4)
Top().mark(bnd,5)
Bottom().mark(bnd,6)
Interface().mark(bnd,7)


# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)
W = FunctionSpace(mesh,'CG', 2)
VQW = MixedFunctionSpace([V,Q,V])

u,p,d = TrialFunctions(VQW)
v,q,w = TestFunctions(VQW)

u_in = Expression(('0.01*(std::exp(-pow((t-0.2),2)/2.0) + std::exp(-pow((t-1.2),2)/2.0) + std::exp(-pow((t-2.2),2)/2.0))','0'),t=0)


p_in = Expression('(std::exp(-pow((t-0.05),2)/0.05) + std::exp(-pow((t-1.05),2)/0.05) + std::exp(-pow((t-2.05),2)/0.05))',t=0)
noslip = Constant((0.0,0.0))
bcd1 = DirichletBC(VQW.sub(2),noslip,bnd,1)
bcd2 = DirichletBC(VQW.sub(2),noslip,bnd,2)
bcd3 = DirichletBC(VQW.sub(2),noslip,bnd,3)
bcd4 = DirichletBC(VQW.sub(2),noslip,bnd,4)
bcd6 = DirichletBC(VQW.sub(2),noslip,bnd,6)

bcu5 = DirichletBC(VQW.sub(0),noslip,bnd,5)
bcu6 = DirichletBC(VQW.sub(0),noslip,bnd,6)
bcu1 = DirichletBC(VQW.sub(0),noslip,bnd,1)
bcu3 = DirichletBC(VQW.sub(0),u_in,bnd,3)
bcu2 = DirichletBC(VQW.sub(0),noslip,bnd,2)

bcp1 = DirichletBC(VQW.sub(1),0,bnd,1)
bcp2 = DirichletBC(VQW.sub(1),0,bnd,2)
bcp3 = DirichletBC(VQW.sub(1),p_in,bnd,	3)
bcp4 = DirichletBC(VQW.sub(1),0,bnd,4)
bcp5 = DirichletBC(VQW.sub(1),0,bnd,5)
bcp6 = DirichletBC(VQW.sub(1),0,bnd,6)


bcd = [bcd1,bcd2,bcd3,bcd4,bcu6]    # skeleton vel
bcu = [bcu3,bcu1,bcu5,bcu6]		# filtration vel
bcp = [bcp2]

bcs = [bu for bu in bcu] + [bd for bd in bcd]

# d - skeleton vel
# u - filtration vel

ds = Measure('ds')[bnd]
dx = Measure('dx')[SD]

eta = Function(V)
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


kappa = Constant(1.4*10**(-15)*(10**6))



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
	#inner(grad(v),sigma_f(u))*dx(0,subdomain_data=SD)

LFW = rho_f/k*inner(u1,v)*dx(0,subdomain_data=SD)

aFV = inner(grad(d),grad(w))*dx(0,subdomain_data=SD)
LFV = inner(noslip,w)*dx(0,subdomain_data=SD)

aFQ = -inner(div(u),q)*dx(0,subdomain_data=SD)
LFQ = Constant(0)*q*dx(0,subdomain_data=SD)

aF = aFV + aFW + aFQ
LF = LFV + LFW + LFQ

a = aS+aF
L = LS+LF

t = dt
count = 0
while t < T + DOLFIN_EPS:
	u_in.t = t
	p_in.t = t
	b = assemble(L)
	eps = 10
	k_iter = 0
	max_iter = 10
	while eps > 1E-4 and k_iter < max_iter:
		A = assemble(a)
		A.ident_zeros()
		[bc.apply(A,b) for bc in bcs]
		solve(A,UPR.vector(),b,'lu')
		u_,p_,d_ = UPR.split(True)
		eps = errornorm(u_,u0,degree_rise=3)
		k_iter += 1
		print 'k: ',k_iter, 'error: %.3e' %eps
		u0.assign(u_)

	print assemble((float(k)*sigma_dev(d_)*n + sigma_dev(eta)*n)[0]*ds(5))
	print assemble(p_*ds(5))
	#bcu[-1] = DirichletBC(VQW.sub(0),d_,bnd,5) # interface
	ufile << u_
	pfile << p_
	dfile << d_
	tfile << eta
	d1.assign(d_)
	d_.vector()[:] *= float(k)
	eta.vector()[:] += d_.vector()[:]
	mesh.move(d_)
	mesh.bounding_box_tree().build(mesh)
	# Move to next time step
	print 't=%.4f'%t
	t += dt
	u1.assign(u_)

	count += 1



