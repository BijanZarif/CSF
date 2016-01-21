from dolfin import *
set_log_active(False)
import sys

h = 60		# heigth
x0 = 5		# inner wall
x1 = 9		# outer wall
s = 1		# central spinal canal

p0 = Point(-x1,0)
p1 = Point(x1,h)


ufile = File("results_mono/velocity.pvd")
pfile = File("results_mono/pressure.pvd")
dfile = File("results_mono/dU.pvd")
tfile = File("results_mono/U.pvd")
#mesh = Mesh('../meshes/cord_w_csc_mm.xml')
mesh = RectangleMesh(p0,p1,18,30)
dt = 0.005
T = 1.5


eps = 1e-10

class InletS(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],h) and on_bnd and abs(x[0]) < x0 + eps
class OutletS(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],0.0) and on_bnd and abs(x[0]) < x0 + eps

class InletFLeft(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],h) and on_bnd and x[0] < - x0 + eps
class OutletFLeft(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],0.0) and on_bnd and x[0] < - x0 + eps

class InletFRight(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],h) and on_bnd and x[0] > x0 - eps
class OutletFRight(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],0.0) and on_bnd and x[0] > x0 - eps

class Walls(SubDomain):
	def inside(self,x,on_bnd):
		return abs(x[0]) > x1 - eps and on_bnd
class Fluid(SubDomain):
        def inside(self,x,on_bounary):
                return x0-eps < abs(x[0]) < x1+eps
class CSC_bnd(SubDomain):
        def inside(self,x,on_bnd):
                xval = s-eps < abs(x[0]) < s+eps
                yval = h/6.0-eps < abs(x[1]) < h*5/6.0 + eps
                xin = -s - eps < x[0] < s + eps
                yb = h/6.0 - eps < x[1] < h/6.0 + eps
                yt = 5*h/6.0 - eps < x[1] < 5*h/6.0 + eps
                return (xin and (yb or yt)) or (yval and xval)
class Interface(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[0],x0) or near(x[0],-x0)

class CSC(SubDomain):
        def inside(self,x,on_bnd):
                return (abs(x[0]) < s + eps and (h/6.0 - eps < x[1] < 5*h/6.0 + eps))


# DOMAIN
SD = MeshFunction('uint',mesh,mesh.topology().dim())
SD.set_all(1)
Fluid().mark(SD,0)
#CSC().mark(SD,0)

bnd = FacetFunction("size_t",mesh)
InletS().mark(bnd,1)
OutletS().mark(bnd,2)
InletFLeft().mark(bnd,3)
OutletFLeft().mark(bnd,4)
InletFRight().mark(bnd,5)
OutletFRight().mark(bnd,6)
Interface().mark(bnd,20	)
Walls().mark(bnd,8)

# TEST AND TRIALFUNCTIONS
V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)
W = FunctionSpace(mesh,'CG', 2)
VQW = MixedFunctionSpace([V,Q,V])

u,p,d = TrialFunctions(VQW)
v,q,w = TestFunctions(VQW)

u_in = Expression(('0.01*(std::exp(-pow((t-0.2),2)/2.0) + std::exp(-pow((t-1.2),2)/2.0) + std::exp(-pow((t-2.2),2)/2.0))','0'),t=0)


u_in_1 = Expression(('0','0.01*sin(2*pi*t)*(x0-x[1])*(x[1]-x1)'),x0=x0,x1=x1,t=0)
u_in_2 = Expression(('0','0.01*sin(2*pi*t)*(x0-x[1])*(x[1]-x1)'),x0=x0,x1=x1,t=0)


p_in = Expression('5*(std::exp(-pow((t-0.05),2)/0.05) + std::exp(-pow((t-1.05),2)/0.05) + std::exp(-pow((t-2.05),2)/0.05))',t=0)
#p_in = Expression(('0.05'),t=0)

noslip = Constant((0.0,0.0))
bcd1 = DirichletBC(VQW.sub(2),noslip,bnd,1)
bcd2 = DirichletBC(VQW.sub(2),noslip,bnd,2)
bcd3 = DirichletBC(VQW.sub(2),noslip,bnd,3)
bcd4 = DirichletBC(VQW.sub(2),noslip,bnd,4)
bcd5 = DirichletBC(VQW.sub(2),noslip,bnd,5)
bcd7 = DirichletBC(VQW.sub(2),noslip,bnd,7)
bcd6 = DirichletBC(VQW.sub(2),noslip,bnd,6)
bcd8 = DirichletBC(VQW.sub(2),noslip,bnd,8)

bcu5 = DirichletBC(VQW.sub(0),u_in_2,bnd,5)
bcu1 = DirichletBC(VQW.sub(0),noslip,bnd,1)
bcu3 = DirichletBC(VQW.sub(0),u_in_1,bnd,3)
bcu2 = DirichletBC(VQW.sub(0),noslip,bnd,2)
bcu7 = DirichletBC(VQW.sub(0),noslip,bnd,7)
bcu8 = DirichletBC(VQW.sub(0),noslip,bnd,8)

bcp1 = DirichletBC(VQW.sub(1),0,bnd,1)
bcp2 = DirichletBC(VQW.sub(1),0,bnd,2)
bcp3 = DirichletBC(VQW.sub(1),p_in,bnd,	3)
bcp4 = DirichletBC(VQW.sub(1),0,bnd,4)
bcp5 = DirichletBC(VQW.sub(1),0,bnd,5)
bcp6 = DirichletBC(VQW.sub(1),0,bnd,6)


bcd = [bcd1,bcd2,bcd8,bcd3,bcd5,bcd4,bcd6]    # skeleton vel
bcu = [bcu2,bcu8,bcu1]		# filtration vel
bcp = [bcp1,bcp2,bcp5]

bcs = [bu for bu in bcu] + [bd for bd in bcd]# + [bp for bp in bcp]

# d - skeleton vel
# u - filtration vel

ds = Measure('ds')[bnd]
dx = Measure('dx')[SD]
dS = Measure('dS')[bnd]

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


kappa = Constant(1.4*10**(-15)*(10**6))



phi = 0.2    # porosity
K_perm = kappa/mu_f

rho_s = Constant(rho_p - rho_f*phi/(1-phi))


def sigma_dev(d):
	return 2*mu_s*sym(grad(d)) + lamda*tr(sym(grad(d)))*Identity(2)

def sigma_f(u):
	return 2*mu_f*sym(grad(u))



aW = rho_p/k*inner(u,v)*dx(1,subdomain_data=SD) + \
	k*inner(grad(v),sigma_dev(u))*dx(1,subdomain_data=SD)

LW = rho_f/k*inner(u1,v)*dx(1,subdomain_data=SD) - \
	inner(grad(v),sigma_dev(eta))*dx(1,subdomain_data=SD)

AV = inner(u,w)*dx(1,subdomain_data=SD) - inner(d,w)*dx(1,subdomain_data=SD)

aS = aW + AV
LS = LW

# FLUID
kappa2 = 10#Constant(1.4*10**(-15)*(10**9))
phi2 = 1    # porosity
K_perm2 = kappa2/mu_f
rho_p2 = Constant(rho_s*(1-phi2) + rho_f*phi)

aFW = rho_f/k*inner(u,v)*dx(0,subdomain_data=SD) + \
	rho_f*inner(grad(u0)*(u+d*k),v)*dx(0,subdomain_data=SD) - \
	 inner(p,div(v))*dx(0,subdomain_data=SD) + \
	mu_f*inner(grad(v),grad(u))*dx(0,subdomain_data=SD)
LFW = rho_f/k*inner(u1,v)*dx(0,subdomain_data=SD) - \
	inner(p_in*n,v)*ds(3) - \
	inner(p_in*n,v)*ds(5)

aFV = k*inner(grad(d),grad(w))*dx(0,subdomain_data=SD)
LFV = -inner(grad(eta),grad(w))*dx(0,subdomain_data=SD)

aFQ = -inner(div(u),q)*dx(0,subdomain_data=SD)
LFQ = Constant(0)*q*dx(0,subdomain_data=SD)
'''
def sigma_dev2(u):
	return mu_s*sym(grad(u)) + 0.0001*lamda*tr(sym(grad(u)))*Identity(2)


aFW = rho_p2/k*inner(d,w)*dx(0,subdomain_data=SD) + \
	rho_f/k*inner(u,w)*dx(0,subdomain_data=SD) + \
	rho_f*inner(grad(u0)*(u-d),w)*dx(0,subdomain_data=SD) - \
	inner(p,div(w))*dx(0,subdomain_data=SD) + \
	mu_f*inner(grad(u),grad(v))*dx(0,subdomain_data=SD) #+ \
	#k*inner(grad(w),sigma_dev2(d))*dx(0,subdomain_data=SD)
LFW = rho_f/k*inner(d1,w)*dx(0,subdomain_data=SD) + \
	rho_f/k*inner(u1,w)*dx(0,subdomain_data=SD)# - \
	#inner(grad(w),sigma_dev2(eta))*dx(0,subdomain_data=SD)

aFV = rho_f/k*inner(d,v)*dx(0,subdomain_data=SD) + \
	(rho_f/(k*phi2)+1./K_perm2)*inner(u,v)*dx(0,subdomain_data=SD) - \
	inner(p,div(v))*dx(0,subdomain_data=SD)
LFV = rho_f/(k*phi2)*inner(u1,v)*dx(0,subdomain_data=SD) + \
	rho_f/k*inner(d1,v)*dx(0,subdomain_data=SD)

aFQ = -inner(div(u+d),q)*dx(0,subdomain_data=SD)
LFQ = Constant(0)*q*dx(0,subdomain_data=SD)
'''
aF = aFV + aFW + aFQ
LF = LFV + LFW + LFQ

a = aS+aF
L = LS+LF

t = dt
count = 0
#plot(bnd)
#interactive()
#sys.exit()

while t < T + DOLFIN_EPS:
	u_in_1.t = t
	u_in_2.t = t
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



