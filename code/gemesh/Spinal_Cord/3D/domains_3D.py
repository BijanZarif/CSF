from dolfin import *
import numpy as np

#mesh = Mesh('meshes/drosdal2D.xml')

#mesh = refine(mesh)
#mesh = refine(mesh)	

rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)


## from Karen/Nina
E = (20*10**3) # Pa --- Ozawa rabbit spinal cord 5*10**3 (No Pia) 16*10^3 (w Pia)
Pr = 0.479#0.479       #

rho_s = Constant(2*rho_f)
lamda = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))

K_perm = Constant(1.4*10**(-15)*(10**6))  # mm^2 use (1.3*10-16 - 1.4*10-14) m^2

#mu_s, lamda_s = mu_s*1000, lamda_s*1000


h = 60		# heigth
x0 = 5		# inner wall
x1 = 9		# outer wall
s = 3		# central spinal canal

tol = 1e-1
tol = 1


class Fluid_out(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[2],0.0) and on_boundary and sqrt(x[0]**2+x[1]**2) > x0 - tol

class Solid_out(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[2],0.0) and on_boundary and sqrt(x[0]**2+x[1]**2) < x0 + tol

class Fluid_walls(SubDomain):
	def inside(self,x,on_boundary):
		return sqrt(x[0]**2+x[1]**2) > x1 - tol and on_boundary

class Interface(SubDomain):
	def inside(self,x,on_boundary):
		return sqrt(x[0]**2+x[1]**2) < x0 + tol and sqrt(x[0]**2+x[1]**2) > x0 - tol

class Fluid_in(SubDomain):
	def inside(self,x,on_boundary):
		return x[2]>(h-tol) and on_boundary and sqrt(x[0]**2+x[1]**2) > x0 - tol

class Solid_in(SubDomain):
	def inside(self,x,on_boundary):
		return x[2]>(h-tol) and on_boundary and sqrt(x[0]**2+x[1]**2) < x0 + tol

class CSC_bnd(SubDomain):
	def inside(self,x,on_boundary):
		xval = s-tol < sqrt(x[0]**2+x[1]**2) < s+tol
		yval = h/6.0-tol < abs(x[2]) < h*5/6.0 + tol
		xin = sqrt(x[0]**2+x[1]**2) < s + tol
		yb = h/6.0 - tol < x[2] < h/6.0 + tol
		yt = 5*h/6.0 - tol < x[2] < 5*h/6.0 + tol
		return (xin and (yb or yt)) or (yval and xval)

class CSC(SubDomain):
	def inside(self,x,on_boundary):
		return (sqrt(x[0]**2+x[1]**2) < s + tol and (h/6.0 - tol < x[2] < 5*h/6.0 + tol))

class Solid(SubDomain):
	def inside(self,x,on_boundary):
		csc = (sqrt(x[0]**2+x[1]**2) < s + tol and (h/6.0 - tol < x[2] < 5*h/6.0 + tol))
		xval = s-tol < sqrt(x[0]**2+x[1]**2) < s+tol
		yval = h/6.0-tol < abs(x[2]) < h*5/6.0 + tol
		xin = sqrt(x[0]**2+x[1]**2) < s + tol
		yb = h/6.0 - tol < x[2] < h/6.0 + tol
		yt = 5*h/6.0 - tol < x[2] < 5*h/6.0 + tol
		csc_bnd = (xin and (yb or yt)) or (yval and xval)
		return sqrt(x[0]**2+x[1]**2) < x0 + tol# and (not csc or csc_bnd)   

class Fluid(SubDomain):
	def inside(self,x,on_bounary):
		return x0-tol < sqrt(x[0]**2+x[1]**2) < x1+tol

class MyEx(Expression):
	def __init__(self,mesh,tau):
		self.mesh = mesh
		self.tau = tau
	def eval_cell(self,values,x,c):
		tau = self.tau(x)
		cell = Cell(self.mesh,c.index)
		n = cell.normal(c.local_facet)
		values[0] = (tau[0]*n[0] + tau[1]*n[1])
		values[1] = (tau[2]*n[0] + tau[3]*n[1])
	def value_shape(self):
		return (2,)


parameters['allow_extrapolation'] = True


