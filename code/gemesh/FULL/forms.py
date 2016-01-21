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

dt = 0.1
T = 1000*dt

mesh = RectangleMesh(p0,p1,L/2,h)

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

plot(SD,wireframe=True)
interactive()

bnd = FacetFunction("size_t",mesh)
InletS().mark(bnd,1)
OutletS().mark(bnd,2)
InletF().mark(bnd,3)
OutletF().mark(bnd,4)
Top().mark(bnd,5)
Bottom().mark(bnd,6)
Interface().mark(bnd,7)

# CONSTANTS
nu = 0.658
rho_f = 1.0*1e-3
mu_f = rho_f*nu
E = 5*10**6
Pr = 0.479
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))
rho_s = Constant(500*rho_f)
kappa = Constant(1.4*10**(-15)*(10**6))
phi = 0.2    									# porosity
K_perm = kappa/mu_f
rho_p = Constant(rho_s*(1-phi) + rho_f*phi)

# SOLID

def LHS_S_mom(u,p,d,v):
	return rho_p/k*inner(u,v)*dx + \
	rho_f/k*inner(d,v)*dx - \
	inner(p,div(v))*dx + k*inner(grad(v),sigma_dev(u))*dx

def RHS_S_mom(u1,d1,v,eta):
	return rho_f/k*inner(u1,v)*dx + \
	rho_f/k*inner(d1,v)*dx - \
	inner(grad(v),sigma_dev(eta))*dx 

def LHS_S_vel(u,p,d,w)
	return rho_f/k*inner(u,w)*dx + \
	(rho_f/(k*phi)+1./K_perm)*inner(d,w)*dx - \
	inner(p,div(w))*dx

def RHS_S_vel(u1,d1):
	return rho_f/(k*phi)*inner(d1,w)*dx + \
	rho_f/k*inner(u1,w)*dx #- \
#	 inner(p_in*n,w)*ds(1) + \
#	inner(p_in*n,w)*ds(2)

def LHS_S_con(u,d,q):
	return -inner(div(u)+div(d),q)*dx

def RHS_S_con(q):
	return Constant(0)*q*dx

# FLUID
