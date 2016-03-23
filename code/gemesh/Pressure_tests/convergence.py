from dolfin import *
from pylab import find,array
set_log_active(False)
import sys


mesh = UnitSquareMesh(20,20)

l = 1
h = 1

w_e = Expression(('0.1*sin(pi*x[0])*sin(pi*t)','-0.1*sin(pi*x[1])*sin(pi*t)'), t=0)
u_e = Expression(('0.1*sin(pi*x[0])*sin(pi*t)','-0.1*sin(pi*x[1])*sin(pi*t)'), t=0)
p_e = Expression(('0.1*sin(pi*x[1])*sin(pi*x[0])*sin(pi*t)'), t=0)


Q = VectorFunctionSpace(mesh,'CG',4)
def epsilon(u):
	return sym(grad(u))

u_t = Expression(('sin(pi*x[0])*cos(t)','-sin(pi*x[1])*cos(t)'), t=0)
dP = Expression(('cos(pi*x[0])*sin(pi*x[1])*sin(t)','cos(pi*x[1])*sin(pi*x[0])*sin(t)'), t=0)

V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
W = VectorFunctionSpace(mesh,'CG',1)
VPW = MixedFunctionSpace([V,P,W])

U = Function(W)

f_f = u_t + dP - 2*div(epsilon(interpolate(u_e,Q)))
f_s = u_t - 2*div(epsilon(interpolate(U,Q)))

v,p,w = TrialFunctions(VPW)
phi,eta,psi = TestFunctions(VPW)

eps = 1e-8

class Walls(SubDomain):
	def inside(self,x,on_bnd):
		return (x[1] > h - eps or x[1	] < eps) and on_bnd
class Inlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[0] > l - eps and on_bnd
class Outlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[0] < eps and on_bnd

class Solid(SubDomain):
	def inside(self,x,on_bnd):
		return x[0] > l/2.0 - eps

class Interface(SubDomain):
	def inside(self,x,on_bnd):
		return abs(x[0] - l/2.) < eps

def boundary(x,on_bnd):
	return on_bnd


bnd = FacetFunction('size_t',mesh)
SD = MeshFunction('size_t',mesh,2)
bnd.set_all(0)
SD.set_all(0)
Walls().mark(bnd,1)
Outlet().mark(bnd,2)
Inlet().mark(bnd,3)
Interface().mark(bnd,4)
Solid().mark(SD,1)
noslip = Constant((0,0))
ds = Measure('ds')[bnd]
dx = Measure('dx')[SD]
dS = Measure('dS')[bnd]
dx_f = dx(0,subdomain_data=SD)
dx_s = dx(1,subdomain_data=SD)

bcs = [DirichletBC(VPW.sub(0),u_e,boundary),DirichletBC(VPW.sub(1),p_e,boundary),DirichletBC(VPW.sub(2),w_e,boundary)]

vfile = File('results_channel/v.pvd')
pfile = File('results_channel/p.pvd')
wfile = File('results_channel/w.pvd')


v0 = Function(V)
v1 = Function(V)
dt = 0.02
k = Constant(dt)


rho_f = Constant(1.)		# g/mm
nu_f = Constant(1)			# mm**2/s
mu_f = Constant(1)		# g/(mm*s)
rho_s = Constant(1)
mu_s = Constant(1)
lamda = Constant(1)
def sigma_dev(U):
	return 2*mu_s*sym(grad(U)) + lamda*tr(sym(grad(U)))*Identity(2)

delta = 1e8
n = FacetNormal(mesh)

aMS = rho_s/k*inner(v,phi)*dx_s + \
	k*inner(sigma_dev(v),grad(phi))*dx_s# - \
#	2*mu_f*inner(-p('+')*Identity(2)*n('+') + epsilon(v('+'))*n('+'),phi('+'))*dS(4)

LMS = rho_s/k*inner(v1,phi)*dx_s - \
	inner(sigma_dev(U),grad(phi))*dx_s + \
	inner(f_s,phi)*dx_s

aDS = delta*inner(v,psi)*dx_s - delta*inner(w,psi)*dx_s


aCS = -inner(div(v),eta)*dx_s


aS = aMS + aDS# + aCS
LS = LMS

# FLUID

aMF = rho_f/k*inner(v,phi)*dx_f + \
	rho_f*inner(grad(v0)*(v-w),phi)*dx_f - \
	 inner(p,div(phi))*dx_f + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx_f# - \
#	2*k*mu_s*inner(sigma_dev(w('-'))*n('-'),phi('-'))*dS(4)

LMF = rho_f/k*inner(v1,phi)*dx_f + \
	  inner(f_f,phi)*dx_f# + \
#	  2*mu_s*inner(sigma_dev(U('-'))*n('-'),phi('-'))*dS(4)

aDF = k*inner(grad(w),grad(psi))*dx_f
LDF = -inner(grad(U),grad(psi))*dx_f

aCF = -inner(div(v),eta)*dx_f
LFQ = Constant(0)*eta*dx_f

aF = aMF + aDF + aCF
LF = LMF + LDF + LFQ

a = aF + aS
L = LF + LS



t = 0
T = 2
dt = 0.01

VPW_ = Function(VPW)


while t < T + DOLFIN_EPS:
	u_e.t = t
	p_e.t = t
	w_e.t = t
	b = assemble(L)
	err = 10
	k_iter = 0
	max_iter = 5
	while err > 1E-6 and k_iter < max_iter:
		A = assemble(a)
		[bc.apply(A,b) for bc in bcs]
		solve(A,VPW_.vector(),b,'lu')
		v_,p_,w_ = VPW_.split(True)
		err = errornorm(v_,v0,degree_rise=3)
		k_iter += 1
		print 'k: ',k_iter, 'error: %.3e' %err
		v0.assign(v_)
	vfile << v_
	pfile << p_	
	wfile << w_	
	w_.vector()[:] *= float(k)
	U.vector()[:] += w_.vector()[:]
	mesh.move(w_)
	mesh.bounding_box_tree().build(mesh)
	v1.assign(v_)

	
	print 't=%.4f'%t
	t += dt
'''
A = assemble(a)
b = assemble(L)
up = Function(VP)
[bc.apply(A,b) for bc in bcs]
solve(A,up.vector(),b)
u_,p_ = up.split(True)
plot(p_)
interactive()
ufile << u_
pfile << p_
'''
