from FSI_Solver import *
set_log_active(False)

class Biot_Solver(FSI_Solver):
	def __init__(self, rho_f, nu_f, E, Poisson, K_perm=Constant(1.4*10**(-15)*(10**6))):
		FSI_Solver.__init__(self,rho_f, nu_f, E, Poisson)
		self.kappa = Constant(K_perm/float(self.nu_f*self.rho_f))
		
	
	def make_functions(self, mesh_file='meshes/wide_syrinx.xml',refine=0):
		
		[[v,w,p],[phi,eta,psi],[dx_f,dx_s]] = FSI_Solver.make_functions(self, mesh_file='meshes/wide_syrinx.xml',refine=0)
		self.a = self.a - inner(p,div(phi))*dx_f-inner(div(v),eta)*dx_s - self.kappa*inner(grad(p),grad(eta))*dx_s


