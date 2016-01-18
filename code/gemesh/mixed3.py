from domains_and_boundaries import *
#mesh = refine(mesh)
#mesh =refine(mesh)

subdomains = MeshFunction("uint", mesh, mesh.topology().dim())
subdomains.set_all(0)
CSF_mesh = SubMesh(mesh,CSF_s)
Cord_mesh = SubMesh(mesh,cord)
CSC_mesh = SubMesh(mesh,in_csc)

CSF_s.mark(subdomains,0)
cord.mark(subdomains,1)
in_csc.mark(subdomains,1)

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

bcs = [bcu1,bcw1,bcw2,bcu2,bcu5,bcp4]

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
D1 = inner(grad(u),grad(r))*dx(1,subdomain_data=subdomains) + \
     inner(f,r)*dx(1,subdomain_data=subdomains) + \
	 inner(r('+'),dot(grad(u('+')),n('+')))*dS(2) - \
	 inner(p0*n,r)*ds(6)

### Biot (momentum + continuity)

B_m = mu_s*inner(grad(w),grad(v))*dx(1,subdomain_data=subdomains) + \
	  lamda_s*inner(div(w),div(v))*dx(1,subdomain_data=subdomains) - \
	  inner(div(v),p)*dx(1,subdomain_data=subdomains) + \
	  inner(p0*n,v)*ds(6) - \
	  inner(p0*n,v)*ds(5)

B_c = -kappa*inner(grad(q),grad(p))*dx(1,subdomain_data=subdomains)


B_t =  1./k*inner(w-w1,grad(q))*dx(1,subdomain_data=subdomains)


F_Cord = D1 + B_m + B_c

up = Constant((1.0,0.0))
down = Constant((-1.0,0.0))

# Dummy for w:
D2 = inner(grad(w),grad(v))*dx(0,subdomain_data=subdomains) + \
     inner(dot(grad(w('-')),n('-')),v('-'))*dS(2) - \
	 inner(up,v)*ds(3) - \
	 inner(down,v)*ds(4)

S_m = mu_f*inner(grad(u),grad(v))*dx(0, subdomain_data=subdomains) - \
      inner(div(v),p)*dx(0, subdomain_data=subdomains) + \
      inner(p0*n,v)*ds(3) - \
      inner(p0*n,v)*ds(4)

S_c = inner(grad(q),u)*dx(0, subdomain_data=subdomains) - \
      inner(q,dot(u,n))*ds(3) - \
	  inner(q,dot(u,n))*ds(4)


F_STOKES = D2 + S_m + S_c


F = F_STOKES + F_Cord

UPW_ = Function(VQR)


solve(lhs(F)==rhs(F), UPW_,bcs=bcs,solver_parameters={"linear_solver": "lu"})
print 'hi'
u,p,w = UPW_.split(True)
for x in w.vector().array():
	print x
plot(subdomains)
interactive()

