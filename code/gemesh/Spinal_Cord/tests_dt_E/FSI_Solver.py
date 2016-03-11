from dolfin import *
#from pylab import array, append, find

h = 60+5		# heigth   #Warning
x0 = 5		# inner wall
x1 = 9		# outer wall
s = 3		# central spinal canal

tol = 1e-12

class Fluid_out(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[1],0.0) and on_boundary and abs(x[0]) > x0 - tol

class Solid_out(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[1],0.0) and on_boundary and abs(x[0]) < x0 + tol

class Fluid_walls(SubDomain):
        def inside(self,x,on_boundary):
                return abs(x[0]) > x1 - tol and on_boundary 

class Interface(SubDomain):
        def inside(self,x,on_boundary):
                return abs(x[0]) < x0 + tol and abs(x[0]) > x0 - tol

class Fluid_in_l(SubDomain):
        def inside(self,x,on_boundary):
                return x[1]>(h-tol) and on_boundary and x[0] < - x0 + tol

class Fluid_in_r(SubDomain):
		def inside(self,x,on_boundary):
				return x[1]>(h-tol) and on_boundary and x[0] > x0 - tol

class Solid_in(SubDomain):
        def inside(self,x,on_boundary):
                return x[1]>(h-tol) and on_boundary and abs(x[0]) < x0 + tol

h = 60 # Warning!

class CSC_bnd(SubDomain):
        def inside(self,x,on_boundary):
                xval = s-tol < abs(x[0]) < s+tol
                yval = h/6.0-tol < abs(x[1]) < h*5/6.0 + tol
                xin = -s - tol < x[0] < s + tol
                yb = h/6.0 - tol < x[1] < h/6.0 + tol
                yt = 5*h/6.0 - tol < x[1] < 5*h/6.0 + tol
                return (xin and (yb or yt)) or (yval and xval)

class CSC(SubDomain):
        def inside(self,x,on_boundary):
                return (abs(x[0]) < s + tol and (h/6.0 - tol < x[1] < 5*h/6.0 + tol))

class Solid(SubDomain):
        def inside(self,x,on_boundary):
                csc = (abs(x[0]) < s + tol and (h/6.0 - tol < x[1] < 5*h/6.0 + tol))
		xval = s-tol < abs(x[0]) < s+tol
                yval = h/6.0-tol < abs(x[1]) < h*5/6.0 + tol
                xin = -s - tol < x[0] < s + tol
                yb = h/6.0 - tol < x[1] < h/6.0 + tol
                yt = 5*h/6.0 - tol < x[1] < 5*h/6.0 + tol
                csc_bnd = (xin and (yb or yt)) or (yval and xval)
                return abs(x[0]) < x0 + tol and (not csc or csc_bnd)   

class Fluid(SubDomain):
        def inside(self,x,on_bounary):
                return x0-tol < abs(x[0]) < x1+tol



class FSI_Solver:
	def __init__(self, rho_f, nu_f, E, Poisson):
		self.rho_f = Constant(rho_f)
		self.nu_f = Constant(nu_f)

		self.rho_s = Constant(2*rho_f)
		self.lamda = Constant(E*Poisson/((1.0+Poisson)*(1.0-2*Poisson)))
		self.mu_s = Constant(E/(2*(1.0+Poisson)))
		
	'''
	def Brucker_inlet(self):
		C1_Brucker = 10*array([2.26, 1.53, 0.18, 0.19, 0.2, 0.05, -0.39, -0.69, -0.93, -0.66, -0.63, -0.7, -0.23, 2.51])      # 12:00

		C1_Brucker2 = 10*array([-0.33, -0.35, -0.51, 0.99, 1.27, 0.83, 0.71, 0.67, -0.15, -0.71, -0.05, -0.21, -0.43, -0.62]) # 6:00


		t_Brucker = 1e-3*array([10, 73, 136, 199, 262, 325, 388, 451, 514, 577, 640, 703, 766, 829])
		C = C1_Brucker[5:]
		C2 = C1_Brucker[:5]
		C1_Brucker = append(C,C2)

		Cop = C1_Brucker2[5:]
		Cop2 = C1_Brucker2[:5]
		C1_Brucker2 = append(Cop,Cop2)

		class MyExpression0(Expression):
			def __init__(self,t_Brucker,C1_Brucker,t):
				self.t = t
				self.t_Brucker = t_Brucker
				self.C1_Brucker = C1_Brucker

			def eval(self,values,x):
				t = self.t
				t_Brucker = self.t_Brucker
				C1_Brucker = self.C1_Brucker
				while t > 0.829: 
					t -= 0.829
				tval = t_Brucker
				yval = C1_Brucker
				idx = find(t_Brucker >= t)[0] - 1

				values[1] = -(yval[idx] + (yval[idx+1]-yval[idx])*(t-tval[idx])/(tval[idx+1]-tval[idx]))
				values[0] = 0

			def value_shape(self):
				return (2,)
		return MyExpression0(t_Brucker,C1_Brucker,0.0)
	'''
	def make_functions(self, mesh_file='meshes/wide_syrinx.xml',refine=0):
		self.mesh = Mesh(mesh_file)
		refinements = 0
		while refine > refinements:
			self.mesh = refine(self.mesh)
			refinements += 1
		
		rho_f, nu_f, rho_s, lamda, mu_s = self.rho_f, self.nu_f, self.rho_s, self.lamda, self.mu_s
		mu_f = rho_f*nu_f
	

		SD = MeshFunction('size_t', self.mesh, self.mesh.topology().dim())
		SD.set_all(0)
		Solid().mark(SD,1)
		#CSC().mark(SD,1)
		
		# DEFINING BOUNDARIES
		boundaries = FacetFunction("size_t",self.mesh)
		boundaries.set_all(0)
		Fluid_in_l().mark(boundaries,1)
		Fluid_in_r().mark(boundaries,2)
		Solid_in().mark(boundaries,3)
		Fluid_out().mark(boundaries,4)
		Solid_out().mark(boundaries,5)
		Interface().mark(boundaries,6)
		Fluid_walls().mark(boundaries,7)
		CSC_bnd().mark(boundaries,8)
		# TEST AND TRIALFUNCTIONS
		V = VectorFunctionSpace(self.mesh,'CG',2)
		P = FunctionSpace(self.mesh,'CG',1)
		self.V = V
		VPW = MixedFunctionSpace([V,P,V])
		self.VPW = VPW

		v,p,w = TrialFunctions(VPW)
		phi,eta,psi = TestFunctions(VPW)


		self.v0 = Function(V)
		self.v1 = Function(V)
		self.U = Function(V)
		
		v0,v1,U = self.v0, self.v1, self.U
		
		# INITIAL AND BOUNDARY CONDITIONS

		# FLUID
		noslip = Constant((0.0,0.0))
		self.vf_l = Constant((0,-1))#self.Brucker_inlet()
		self.vf_r = Constant((0,-1))#self.Brucker_inlet()
		bcv1 = DirichletBC(VPW.sub(0),self.vf_l,boundaries,1)     # Fluid in_l
		bcv2 = DirichletBC(VPW.sub(0),self.vf_r,boundaries,2)	# Fluid in_r
		bcv3 = DirichletBC(VPW.sub(0),noslip,boundaries,3) # Solid in
		bcv4 = DirichletBC(VPW.sub(0),noslip,boundaries,4) # Fluid out
		bcv5 = DirichletBC(VPW.sub(0),noslip,boundaries,5) # Solid out
		bcv6 = DirichletBC(VPW.sub(0),noslip,boundaries,6) # Interface
		bcv7 = DirichletBC(VPW.sub(0),noslip,boundaries,7) # Fluid walls


		bcv = [bcv1, bcv2, bcv3, bcv5, bcv7] # don't use bcv6 for FSI

		# SOLID

		# MESH DISPLACEMENT

		bcw1 = DirichletBC(VPW.sub(2),noslip,boundaries,1)  # Fluid in_l
		bcw2 = DirichletBC(VPW.sub(2),noslip,boundaries,2)  # Fluid in_r
		bcw3 = DirichletBC(VPW.sub(2),noslip,boundaries,3)  # Solid in
		bcw4 = DirichletBC(VPW.sub(2),noslip,boundaries,4)  # Fluid out
		bcw5 = DirichletBC(VPW.sub(2),noslip,boundaries,5)  # Solid out
		bcw6 = DirichletBC(VPW.sub(2),noslip,boundaries,6)  # Interface
		bcw7 = DirichletBC(VPW.sub(2),noslip,boundaries,7) # Fluid walls
		bcw = [bcw1,bcw2,bcw3,bcw4,bcw5,bcw7]

		# CREATE FUNCTIONS
		self.bcs = bcv + bcw

		# Define coefficients
		n = FacetNormal(self.mesh)


		dS = Measure('dS')[boundaries]
		dx = Measure('dx')[SD]
		ds = Measure('ds')[boundaries]

		dx_f = dx(0,subdomain_data=SD)
		dx_s = dx(1,subdomain_data=SD)


		def sigma_dev(U):
			return 2*mu_s*sym(grad(U)) + lamda*tr(sym(grad(U)))*Identity(2)

		def sigma_f(v):
			return 2*mu_f*sym(grad(v))
	
		self.dt = 0.01
		k = self.dt

		epsilon = 1e8

		aMS = rho_s/k*inner(v,phi)*dx_s + \
			k*inner(sigma_dev(v),grad(phi))*dx_s

		LMS = rho_s/k*inner(self.v1,phi)*dx_s - \
			inner(sigma_dev(U),grad(phi))*dx_s# - \
			#inner(g,phi)*dx_s

		aDS = epsilon*inner(v,psi)*dx_s - epsilon*inner(w,psi)*dx_s


		aCS = -inner(div(v),eta)*dx_s - epsilon*inner(p,eta)*dx_s


		aS = aMS + aDS# + aCS
		LS = LMS

		# FLUID
		k = Constant(self.dt)

		aMF = rho_f/k*inner(v,phi)*dx_f + \
			rho_f*inner(grad(self.v0)*(v-w),phi)*dx_f - \
			 inner(p,div(phi))*dx_f + \
			2*mu_f*inner(sym(grad(v)),sym(grad(phi)))*dx_f

		LMF = rho_f/k*inner(self.v1,phi)*dx_f

		aDF = k*inner(grad(w),grad(psi))*dx_f
		LDF = -inner(grad(self.U),grad(psi))*dx_f

		aCF = -inner(div(v),eta)*dx_f
		LFQ = Constant(0)*eta*dx_f

		aF = aMF + aDF + aCF
		LF = LMF + LDF + LFQ

		a = aS+aF
		L = LS+LF
		self.a = a
		self.L = L

	def run(self,dt,T,folder,save_count=5,solver=None, prec=None,max_iter=6,max_error=1E-6):
		ufile = File("%s/velocity.xdmf"%folder) # xdmf
		pfile = File("%s/pressure.xdmf"%folder)
		dfile = File("%s/dU.xdmf"%folder)
		tfile = File("%s/U.xdmf"%folder)
		VPW_ = Function(self.VPW)
		self.v0.vector().array()[:] = 0
		self.v1.vector().array()[:] = 0
		self.dt = dt                 # To update k in a, L
		t = dt
		count = 0
		while t < T + DOLFIN_EPS:# and (abs(FdC) > 1e-3 or abs(FlC) > 1e-3):
				self.vf_l.t=t
				self.vf_r.t=t
				b = assemble(self.L)
				eps = 10
				k_iter = 0
				while eps > max_error and k_iter < max_iter:
					A = assemble(self.a)
					A.ident_zeros()
					[bc.apply(A,b) for bc in self.bcs]
					if solver is None:
						solve(A,VPW_.vector(),b)
						print 'LU'
					elif prec is None:
						solve(A,VPW_.vector(),b, solver)
					else:
						solve(A,VPW_.vector(),b, solver, prec)
					v_,p_,w_ = VPW_.split(True)
					eps = errornorm(v_,self.v0,degree_rise=3)
					k_iter += 1
					print 'k: ',k_iter, 'error: %.3e' %eps
					self.v0.assign(v_)

				if count%save_count==0:
					ufile << v_
					pfile << p_
					dfile << w_
					tfile << self.U
				w_.vector()[:] *= dt
				self.U.vector()[:] += w_.vector()[:]
				self.mesh.move(w_)
				self.mesh.bounding_box_tree().build(self.mesh)
				self.v1.assign(v_)
				print 't=%.4f'%t
				t += dt
				count += 1


