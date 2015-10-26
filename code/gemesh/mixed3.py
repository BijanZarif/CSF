from domains_and_boundaries import *
#mesh = refine(mesh)
#mesh =refine(mesh)

subdomains = MeshFunction("uint", mesh, mesh.topology().dim())

CSF_mesh = SubMesh(mesh,CSF_s)
Cord_mesh = SubMesh(mesh,cord)
CSC_mesh = SubMesh(mesh,in_csc)

CSF_s.mark(subdomains,0)
cord.mark(subdomains,1)
#in_csc.mark(subdomains,0)

#mesh = CSF_mesh

boundaries = FacetFunction('size_t', mesh)#MeshFunction('uint',mesh ,mesh.topology().dim()-1)
boundaries.set_all(0)

out_wall.mark(boundaries,1)
in_wall.mark(boundaries,2)
top_out.mark(boundaries,3)
bot_out.mark(boundaries,4)
bot_in.mark(boundaries,5)
top_in.mark(boundaries,6)
csc_bnd.mark(boundaries,7)

"""
plot(mesh)
interactive()
import sys
sys.exit()
"""

V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)
W = TensorFunctionSpace(mesh,'CG',2)
VQR = MixedFunctionSpace([V,Q,V])
upw = TrialFunction(VQR)
v,q,r = TestFunctions(VQR)


noslip = Constant((0.0,0.0))

pb = Constant(0.0)
P = 1
pt = Expression(('P*(cos(2*pi*t) + 0.2*cos(8*pi*t))'),t=0,P=P)

flow = Expression(('0.0','P*(cos(2*pi*t) + 0.2*cos(8*pi*t))*(std::abs(x[0])-x0)*(std::abs(x[0])-x1)'),x0=x0,x1=x1,t=0,P=P)

p0 = pt#interpolate(pt,Q)


bcu1 = DirichletBC(VQR.sub(0),noslip,boundaries, 1)  # outer wall
bcu2 = DirichletBC(VQR.sub(0),noslip,boundaries, 2)  # inner wall
bcu5 = DirichletBC(VQR.sub(0),noslip,boundaries, 5)  # bottom cord
bcu6 = DirichletBC(VQR.sub(0),noslip,boundaries, 6)  # top cord
bcu7 = DirichletBC(VQR.sub(0),noslip,boundaries, 7)  # CSC boundary

bcw1 = DirichletBC(VQR.sub(2),noslip,boundaries,1)
bcw2 = DirichletBC(VQR.sub(2),noslip,boundaries,2)  # inner wall
bcw7 = DirichletBC(VQR.sub(2),noslip,boundaries,7)  # CSC boundary
bcw5 = DirichletBC(VQR.sub(2),noslip,boundaries,5)  # bottom cord
bcw6 = DirichletBC(VQR.sub(2),noslip,boundaries,6)  # top cord
bcw3 = DirichletBC(VQR.sub(2),noslip,boundaries,3)  # top CSF
bcw4 = DirichletBC(VQR.sub(2),noslip,boundaries,4)  # bottom CSF

inflow = DirichletBC(VQR.sub(0),flow,boundaries, 3)  # top CSF

bcU0 = DirichletBC(VQR.sub(0),noslip,boundaries,4)   # bot CSF
bcU1 = DirichletBC(VQR.sub(0),noslip,boundaries,3)   # top CSF

bcp1 = DirichletBC(VQR.sub(1),pt,boundaries, 3)      # top CSF
bcp2 = DirichletBC(VQR.sub(1),pb, boundaries ,4)     # bot CSF
bcp3 = DirichletBC(VQR.sub(1),pt, boundaries, 6)     # top cord
bcp4 = DirichletBC(VQR.sub(1),pb, boundaries, 5)     # bottom cord

dx = Measure("dx")[subdomains]
ds = Measure("ds")[boundaries]
dS = Measure("dS")[boundaries]


bcs = [bcu1,bcw1]

ufile = File('results_mix/vel.pvd')
pfile = File('results_mix/pr.pvd')
wfile = File('results_mix/disp.pvd')

t = 0
T = 0.005
dt = 0.005


# Define coefficients
k = Constant(dt)
f = Constant((0, 0))


n = FacetNormal(mesh)

u,p,w = split(upw)

#up1 = Function(VQR)#,'initial_data/up_mixed.xml')
up_ = Function(V)
u1 = Function(V)#,'initial_data/u_mixed.xml')
p1 = Function(Q)#,'initial_data/p_mixed.xml')
w1 = Function(V)

#sign in third line?
# + inner(q,div(u))*dx(1,subdomain_data=subdomains)
#
kappa = Constant((K_perm/mu_f))

def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)

def sigma(u,p,mu):
    return 2*mu*epsilon(u) - p*Identity(2)


# Dummy for u:
D1 = inner(u,r)*dx(1,subdomain_data=subdomains) + \
     inner(Constant((0,0)),r)*dx(1,subdomain_data=subdomains)

### Biot (momentum + continuity)

B_m = mu_s*inner(grad(w),grad(v))*dx(1,subdomain_data=subdomains) + \
	  lamda_s*inner(div(w),div(v))*dx(1,subdomain_data=subdomains) - \
	  inner(div(v),p)*dx(1,subdomain_data=subdomains) + \
	  inner(p0*n,v)*ds(6) - \
	  inner(p0*n,v)*ds(5)

B_c = kappa*inner(grad(q),grad(p))*dx(1,subdomain_data=subdomains) - \
      inner(q('+'),dot(u('+'),n('+')))*dS(2)


B_t =  1./k*inner(w-w1,grad(q))*dx(1,subdomain_data=subdomains)


F_Cord = D1 + B_m + B_c


# Dummy for w:
D2 = inner(grad(w),grad(v))*dx(0,subdomain_data=subdomains) + \
     inner(dot(grad(w('-')),n('-')),v('-'))*dS(2)

S_m = mu_f*inner(grad(u),grad(v))*dx(0, subdomain_data=subdomains) - \
      inner(div(v),p)*dx(0, subdomain_data=subdomains) + \
      inner(p0*n,v)*ds(3) - \
      inner(p0*n,v)*ds(4)

S_c = inner(grad(q),u)*dx(0, subdomain_data=subdomains) - \
      kappa*inner(q('-'),dot(grad(p('-')),n('-')))*dS(2)# - \
      #1./k*inner(dot(w('-')-w1,n('-')),q('-'))*dS(2)
      


F_STOKES = D2 + S_m + S_c


F = F_STOKES + F_Cord



L = lhs(F)
R = rhs(F)
UPW = Function(VQR)


solve(F == 0, UPW,bcs=bcs,solver_parameters={"linear_solver": "lu"})
print 'hi'
u_,p_,w_ = split(UPW)
plot(UPW)
interactive()
ufile << u_
pfile << p_
wfile << w_
"""
F = action(F,up_)
J = derivative(F, up_, up)
problem = NonlinearVariationalProblem(F,up_,bcs,J)
solver = NonlinearVariationalSolver(problem)

import time

t0 = time.time()


while t<T:
    #pt.t = t
    flow.t = t
    solver.solve()
    u_,p_ = up_.split(True)
    u1.vector()[:] = u_.vector()
    p1.vector()[:] = p_.vector()
    t += dt
    print 't = %g' %t
    ufile << u_
    pfile << p_

V0 = FunctionSpace(mesh, 'DG', 0)
k  = Function(V0)
k_values = [0, 1]
help = np.asarray(subdomains.array(), dtype=np.int32)
k.vector()[:] = np.choose(help, k_values)

U = u_*(k-1) + K*grad(p_)*k
U = project(U,V)

plot(U[0])
interactive()


File('initial_data/u_mixed.xml') << u_
File('initial_data/p_mixed.xml') << p_

print 'time spent: ', time.time()-t0


while t<T:
    pt.t = t

    A = assemble(a)
    b = assemble(L)
    [bc.apply(A,b) for bc in bcs]
    solve(A,up.vector(),b)
    
    up1.assign(up)
    u_,p_ = up.split(True)
    u1.assign(u_)
    t += dt
    print 't = %g' %t
    ufile << u1
    pfile << p_

"""
"""
#a = (1./k)*inner(u,v)*dx + nu_f*inner(grad(u),grad(v))*dx - 1./rho_f*inner(div(v),p)*dx + inner(grad(u)*u1,v)*dx - inner(q,div(u))*dx

#L = rho_f*(1./k)*inner(u1,v)*dx #+ inner(Constant('0'),q)*dx


F_Cord = 2*mu_s*inner(epsilon(u),epsilon(v))*dx(1,subdomain_data=subdomains) + \
         (mu_s+lamda_s)*inner(div(u),div(v))*dx(1,subdomain_data=subdomains) - \
         inner(div(v),p)*dx(1,subdomain_data=subdomains) + \
         inner(p0*n,v)*ds(6) + \
         K*inner(grad(q),grad(p))*dx(1,subdomain_data=subdomains)  + \
         inner(q,dot(u,n))*ds(2) + \
         inner(q,dot(u,n))*ds(7) - \
         K*inner(q,dot(grad(p),n))*ds(5) - \
         K*inner(q,dot(grad(p),n))*ds(6)

F_STOKES = 2*mu_f*inner(epsilon(u),epsilon(v))*dx(0, subdomain_data=subdomains) - inner(div(v),p)*dx(0, subdomain_data=subdomains) - inner(q,div(u))*dx(0, subdomain_data=subdomains)  + inner(p0*n,v)*ds(3) - mu*inner(v,dot(n,grad(u).T))*ds(3) - mu*inner(v,dot(n,grad(u).T))*ds(3)
"""
