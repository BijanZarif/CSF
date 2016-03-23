from dolfin import *
from pylab import find,array, linspace
set_log_active(False)
import scipy.interpolate as Spline

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
		values[0] = self.amp*Spline.UnivariateSpline(tval,yval)(t)

h = 4
l = 60

P1 = Point(0,0)
P2 = Point(l,h)

Nx = 60
Ny = 8

mesh = RectangleMesh(P1,P2,Nx,Ny)

E = Expression(('2*x[0]','0'))

V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
VP = V*P

u,p = TrialFunctions(VP)
v,q = TestFunctions(VP)

eps = 1e-8

class Walls(SubDomain):
	def inside(self,x,on_bnd):
		return (x[1] > h - eps or x[1] < eps) and on_bnd
class Inlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[0] < eps and on_bnd
class Outlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[0] > l - eps and on_bnd


bnd = FacetFunction('size_t',mesh)
Walls().mark(bnd,1)
Outlet().mark(bnd,2)
Inlet().mark(bnd,3)
noslip = Constant((0,0))
ds = Measure('ds')[bnd]
bcs = [DirichletBC(VP.sub(0),noslip,bnd,1),DirichletBC(VP.sub(1),Constant(0),bnd,2),DirichletBC(VP.sub(1),Constant(1),bnd,3)]


ufile = File('results_channel/v.pvd')
pfile = File('results_channel/p.pvd')

VP_ = Function(VP)
u0 = Function(V)
u1 = Function(V)
dt = 0.001
k = Constant(dt)

time1 = array([0, 0.16, 1.1])
press1 = array([2.9, 10, 2.9])
time2 = array([0, 0.3, 1.1])
press2 = array([10, 11.5, 10])
pressure1 = MyExpression0(time1,-press1,dt)
pressure2 = MyExpression0(time2,press2,dt)

t_Erika = linspace(0,1.1,23)
P_Erika = array([-0.011,-0.03,-0.02, 0.002,-0.001, -0.002,-0.003, -0.004,0.001, 0.002,0.003, 0.003, 0.004, 0.004,0.003, 0.004,0.006, 0.04,0.045, 0.01, -0.01,-0.01,-0.01])
P_Erika *= 133*6
pressure = MyExpression0(t_Erika,P_Erika,0.0,amp=1)

pressure2 = Expression(('-10*sin(2*pi*t)'),t=0)
rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)

n = FacetNormal(mesh)
def epsilon(u):
	return sym(grad(u))


a = rho_f*1/k*inner(u,v)*dx + rho_f*inner(grad(u0)*u,v)*dx + 2*mu_f*inner(epsilon(u),epsilon(v))*dx + inner(v,grad(p))*dx - inner(div(u),q)*dx# - mu_f*inner(grad(u).T*n,v)*ds(3) - inner(grad(u).T*n,v)*ds(2)
L = rho_f*1/k*inner(u1,v)*dx# - inner(pressure*n,v)*ds(3)# - inner(pressure2*n,v)*ds(2)

t = dt
T = 2

while t < T + DOLFIN_EPS:
	#if t < 2.0:
	#	pressure1.amp = t/2.
	#pressure.t=t
	b = assemble(L)
	err = 10
	k_iter = 0
	max_iter = 8
	while err > 1E-10 and k_iter < max_iter:
		A = assemble(a)
		[bc.apply(A,b) for bc in bcs]
		solve(A,VP_.vector(),b,'lu')
		u_,p_ = VP_.split(True)
		err = errornorm(u_,u0,degree_rise=3)
		k_iter += 1
		u0.assign(u_)
		print 'k: ',k_iter, 'error: %.3e' %err
	ufile << u_
	pfile << p_
	u1.assign(u_)
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
