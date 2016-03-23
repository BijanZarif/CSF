from dolfin import *
from pylab import find, array
set_log_active(False)
class MyExpression0(Expression):
	def __init__(self,t_Brucker,C1_Brucker,t):
		self.t = t
		self.t_Brucker = t_Brucker
		self.C1_Brucker = C1_Brucker

	def eval(self,values,x):
		t = self.t
		t_Brucker = self.t_Brucker
		C1_Brucker = self.C1_Brucker
		while t > self.t_Brucker[-1]:
			t -= self.t_Brucker[-1]
		tval = t_Brucker
		yval = C1_Brucker
		idx = find(t_Brucker >= t)[0] - 1
		values[1] = -(yval[idx] + (yval[idx+1]-yval[idx])*(t-tval[idx])/(tval[idx+1]-tval[idx]))
		values[0] = 0
	def value_shape(self):
		return (2,)

h = 450
l = 10

P1 = Point(0,0)
P2 = Point(l,h)


mesh = RectangleMesh(P1,P2,5,90)

V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
VP = MixedFunctionSpace([V,P])

v,p = TrialFunctions(VP)
phi,eta = TestFunctions(VP)

eps = 1e-8

class Inlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[1] > h - eps and on_bnd

class Walls(SubDomain):
	def inside(self,x,on_bnd):
		return (x[0] > l - eps or x[0] < eps) and on_bnd
class Outlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[1] < eps and on_bnd


bnd = FacetFunction('size_t',mesh)
Walls().mark(bnd,0)
Inlet().mark(bnd,1)
Outlet().mark(bnd,2)
plot(bnd)
interactive()
noslip = Constant((1,0))

bcs = [DirichletBC(VP.sub(0),noslip,bnd,0)]
ds = Measure('ds')[bnd]

rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)

time1 = array([0, 0.16, 1.1])
press1 = array([2.9, 10, 2.9])
time2 = array([0, 0.3, 1.1])
press2 = array([10, 11.5, 10])

pressure1 = MyExpression0(time1,press1,0)
pressure2 = MyExpression0(time2,press2,0)


dt = 0.01
k = Constant(dt)
v0 = Function(V)
v1 = Function(V)


a = rho_f/k*inner(v,phi)*dx + \
	rho_f*inner(grad(v0)*v,phi)*dx - \
	inner(p,div(phi))*dx - \
	inner(eta,div(v))*dx + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx

L = rho_f/k*inner(v1,phi)*dx# - \
#	  inner(Constant((1,0)),phi)*ds(1) - \
#	  inner(Constant((0,0)),phi)*ds(2)

t = dt
T = 10*dt
VP_ = Function(VP)
ufile = File('results_channel/v.pvd')
pfile = File('results_channel/p.pvd')

while t < T + DOLFIN_EPS:
	pressure1.t=t
	pressure2.t=t
	b = assemble(L)
	err = 10
	k_iter = 0
	max_iter = 5
	while err > 1E-6 and k_iter < max_iter:
	    A = assemble(a)
	    [bc.apply(A,b) for bc in bcs]
	    solve(A,VP_.vector(),b,'lu')
	    v_,p_ = VP_.split(True)
	    err = errornorm(v_,v0,degree_rise=3)
	    k_iter += 1
	    print 'k: ',k_iter, 'error: %.3e' %eps
	    v0.assign(v_)
	ufile << v_
	pfile << p_
	v1.assign(v_)
        print 't=%.4f'%t
	t += dt


