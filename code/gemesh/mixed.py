from domains_and_boundaries import *
#mesh = refine(mesh)
mesh =refine(mesh)

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
#csc_bnd.mark(boundaries,7)

V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',2)
W = TensorFunctionSpace(mesh,'CG',1)
VQ = V*Q
up = TrialFunction(VQ)
v, q = TestFunctions(VQ)


noslip = Constant((0.0,0.0))
P = 1
p0 = Expression(('P*(cos(2*pi*t) + 0.2*cos(8*pi*t))'),t=0,P=P)
pb = Constant(0)
flow = Expression(('0.0','P*(cos(2*pi*t) + 0.2*cos(8*pi*t))*(std::abs(x[0])-x0)*(std::abs(x[0])-x1)'),x0=x0,x1=x1,t=0,P=P)


### BOUNDARY CONDITIONS

bcu0 = DirichletBC(VQ.sub(0),noslip,boundaries, 1)  # outer wall
bcu1 = DirichletBC(VQ.sub(0),noslip,boundaries, 2)  # inner wall
bcu2 = DirichletBC(VQ.sub(0),noslip,boundaries, 5)  # bottom cord
bcu3 = DirichletBC(VQ.sub(0),noslip,boundaries, 6)  # top cord
bcu4 = DirichletBC(VQ.sub(0),noslip,boundaries, 7)  # CSC boundary

inflow = DirichletBC(VQ.sub(0),flow,boundaries, 3)  # top CSF

bcU0 = DirichletBC(VQ.sub(0),noslip,boundaries,4)   # bot CSF
bcU1 = DirichletBC(VQ.sub(0),noslip,boundaries,3)   # top CSF

bcp1 = DirichletBC(VQ.sub(1),p0,boundaries, 3)      # top CSF
bcp2 = DirichletBC(VQ.sub(1),pb, boundaries ,4)     # bot CSF
bcp3 = DirichletBC(VQ.sub(1),p0, boundaries, 6)     # top cord
bcp4 = DirichletBC(VQ.sub(1),pb, boundaries, 5)     # bottom cord


bcs = [bcu0]

###

dx = Measure("dx")[subdomains]
ds = Measure("ds")[boundaries]
dS = Measure("dS")[boundaries]


ufile = File('results_mix/vel.pvd')
pfile = File('results_mix/pr.pvd')
dfile = File('results_mix/disp.pvd')

t = 0
T = 0.005
dt = 0.005


# Define coefficients
k = Constant(dt)
n = FacetNormal(mesh)



u,p = split(up)
up_ = Function(VQ)
u1 = Function(V)#,'initial_data/u_mixed.xml')
p1 = Function(Q)#,'initial_data/p_mixed.xml')


kappa = Constant((K_perm/mu_f))


print float(kappa)
def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)

def sigma(u,p,mu):
    return 2*mu*epsilon(u) - p*Identity(2)

### Biot (momentum + continuity)

B_m = mu_s*inner(grad(u),grad(v))*dx(1,subdomain_data=subdomains) + \
      (lamda_s+mu_s)*inner(div(u),div(v))*dx(1,subdomain_data=subdomains) - \
      inner(div(v),p)*dx(1,subdomain_data=subdomains) + \
      inner(p0*n,v)*ds(6) - \
      inner(p0*n,v)*ds(5)

B_c = kappa*inner(grad(q),grad(p))*dx(1,subdomain_data=subdomains) #+ \
      #inner(q('-'),dot(u('-'),n('-')))*dS(2) + \
      #inner(q('-'),dot(u('-'),n('-')))*dS(7)

B_2 =  inner(q,dot(u,n))*ds(5) + \
       inner(q,dot(u,n))*ds(6)

B_3 = -kappa*inner(q,dot(grad(p),n))*ds(5) - \
      kappa*inner(q,dot(grad(p),n))*ds(6)

B_t = 1./k*inner(dot(u-u1,n),q)*ds(5) + \
      1./k*inner(dot(u-u1,n),q)*ds(6)
      

B_c += B_t


# should be ('-')???

# Time-dependent term(s)

B_t =  1./k*inner(u-u1,grad(q))*dx(1,subdomain_data=subdomains)


F_Cord = B_m + B_c + B_t


### STOKES (momentum + continuity)

S_m = rho_f*inner(grad(u)*u1,v)*dx(0,subdomain_data=subdomains) + \
      mu_f*inner(grad(u),grad(v))*dx(0, subdomain_data=subdomains) - \
      inner(div(v),p)*dx(0, subdomain_data=subdomains) + \
      inner(p0*n,v)*ds(3) - \
      inner(p0*n,v)*ds(4)

S_c1 = inner(div(u),q)*dx(0, subdomain_data=subdomains)

S_c = -inner(u,grad(q))*dx(0,subdomain_data=subdomains) + \
      inner(q,dot(u,n))*ds(3) + \
      inner(q,dot(u,n))*ds(4) #- \
      #kappa*inner(q('+'),dot(grad(p('+')),n('+')))*dS(2) - \
      #kappa*inner(q('+'),dot(grad(p('+')),n('+')))*dS(7)
        



# Time-dependent term(s)

S_t = rho_f/k*inner(u-u1,v)*dx

F_STOKES = S_m + S_c#+ S_t


F = F_STOKES + F_Cord



L = lhs(F)
R = rhs(F)
UP = Function(VQ)

solve(L == R, UP, bcs=bcs)
u_,p_ = UP.split(True)
ufile << u_
pfile << p_


x = Expression('x[0]')
y = Expression('x[1]')
x = project(x,Q)
print assemble(inner(x('+'),x('+'))*dS(2)),assemble(inner(x('-'),x('-'))*dS(2))
val = assemble(inner(x('+'),dot(u_('-'),n('-')))*dS(2))
val2 = assemble(inner(x('+'),dot(kappa*grad(p_('-')),n('-')))*dS(2))
print 'assembled dot(u,n)*x on interface: ',val
print 'assembled dot(kappa*grad(p), n): ',val2
import sys
sys.exit()


'''
d = Function(V)
d = interpolate(u_,V)
u = TestFunction(V)
v = TrialFunction(V)

bc = DirichletBC(V,noslip,boundaries,1)
bc2 = DirichletBC(V,d,boundaries,2)
bc3 = DirichletBC(V,d,boundaries,7)
bc4 = DirichletBC(V,d,boundaries,5)
bc5 = DirichletBC(V,d,boundaries,6)

F = inner(grad(u),grad(v))*dx(0,subdomain_data=subdomains) - \
    inner(Constant((0,0)),v)*dx(0,subdomain_data=subdomains) + \
    inner(grad(u),grad(v))*dx(1,subdomain_data=subdomains) - \
    inner(grad(d),grad(v))*dx(1,subdomain_data=subdomains)

u_ = Function(V)

solve(lhs(F)==rhs(F),u_,bcs=[bc,bc2,bc3,bc4,bc5])

dfile << u_

#plot(mesh)
mesh.move(u_)
mesh.bounding_box_tree().build(mesh)
plot(boundaries)
interactive()
'''
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
