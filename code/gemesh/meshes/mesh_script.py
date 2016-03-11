from dolfin import *
import sys

mesh = Mesh('jvmesh.xml')
V = VectorFunctionSpace(mesh,'CG',2)
P = FunctionSpace(mesh,'CG',1)
ufile = File('v.pvd')

VPW = MixedFunctionSpace([V,P])
print VPW.dim()

v,p = TrialFunctions(VPW)
phi,eta = TestFunctions(VPW)

rho_f = (1./1000)		# g/mm
nu_f = (0.658)			# mm**2/s
mu_f = (nu_f*rho_f)		# g/(mm*s)

E = 62.5*10**3
Pr = 0.479

rho_s = Constant(2*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))

v0 = Function(V)
v1 = Function(V)
U = Function(V)
dt = 0.05
k = Constant(dt)

aMF = rho_f/k*inner(v,phi)*dx + \
	rho_f*inner(grad(v0)*(v),phi)*dx - \
	inner(p,div(phi))*dx + \
	2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx

LMF = rho_f/k*inner(v1,phi)*dx


aCF = -inner(div(v),eta)*dx
LFQ = Constant(0)*eta*dx

aF = aMF + aCF
LF = LMF + LFQ

VPW_ = Function(VPW)

eps = 1e-14

class Outlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[1] < 286.107 and on_bnd
class Inlet(SubDomain):
	def inside(self,x,on_bnd):
		return x[1] > 340.42 and on_bnd
class Boundary(SubDomain):
	def inside(self,x,on_bnd):
		return on_bnd
print 'hi'
mf = FacetFunction('size_t',mesh)
Boundary().mark(mf,2)
Outlet().mark(mf,0)
Inlet().mark(mf,1)

noslip = Constant((0,0,0))
inlet = Constant((0,-1,0))
#noslip = Constant((0,0))
#inlet = Constant((1,0))

bc0 = DirichletBC(VPW.sub(0),noslip,mf,2)
bc1 = DirichletBC(VPW.sub(0),inlet,mf,1)
bcs = [bc1,bc0]

t = dt
T = 20*dt

while t < T:
	A = assemble(aF)
	b = assemble(LF)
	[bc.apply(A,b) for bc in bcs]
	solve(A, VPW_.vector(), b)#, 'cg','amg')
	v_,p_ = VPW_.split(True)
	v0.assign(v1)
	v1.assign(v_)
	ufile << v_	
	t+=dt
	print 't=',t


