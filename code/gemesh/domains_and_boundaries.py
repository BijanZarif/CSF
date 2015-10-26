from dolfin import *
import numpy as np
set_log_active(False)
mesh = Mesh('meshes/cord_w_csc_mm.xml')

mesh = Mesh('meshes/drosdal2D.xml')

#mesh = refine(mesh)
#mesh = refine(mesh)	

rho_f = Constant(1./1000)		# g/mm
nu_f = Constant(0.658)			# mm**2/s
mu_f = Constant(nu_f*rho_f)		# g/(mm*s)


## from Karen/Nina
E = (5.0*10**3) # Pa --- Ozawa rabbit spinal cord 5*10**3 (No Pia) 16*10^3 (w Pia)
Pr = 0.479       #

lamda_s = Constant(E*Pr/((1.0+Pr)*(1.0-2*Pr)))
mu_s = Constant(E/(2*(1.0+Pr)))

K_perm = Constant(1.4*10**(-15)*(10**6))  # mm^2 use (1.3*10-16 - 1.4*10-14) m^2

#mu_s, lamda_s = mu_s*1000, lamda_s*1000

print 'mu: ', float(mu_s)
print 'lambda: ', float(lamda_s)

h = 60		# heigth
x0 = 5		# inner wall
x1 = 9		# outer wall
s = 1		# central spinal canal

tol = 1e-12

class bottom_outer(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[1],0.0) and on_boundary and abs(x[0]) > x0 - tol

class bottom_inner(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[1],0.0) and on_boundary and abs(x[0]) < x0 + tol

class outer_wall(SubDomain):
        def inside(self,x,on_boundary):
                return abs(x[0]) > x1 - tol and on_boundary 

class inner_wall(SubDomain):
        def inside(self,x,on_boundary):
                return abs(x[0]) < x0 + tol and abs(x[0]) > x0 - tol

class top_outer(SubDomain):
        def inside(self,x,on_boundary):
                return x[1]>(h-tol) and on_boundary and abs(x[0]) > x0 - tol

class top_inner(SubDomain):
        def inside(self,x,on_boundary):
                return x[1]>(h-tol) and on_boundary and abs(x[0]) < x0 + tol

class csc_boundary(SubDomain):
        def inside(self,x,on_boundary):
                xval = s-tol < abs(x[0]) < s+tol
                yval = h/6.0-tol < abs(x[1]) < h*5/6.0 + tol
                xin = -s - tol < x[0] < s + tol
                yb = h/6.0 - tol < x[1] < h/6.0 + tol
                yt = 5*h/6.0 - tol < x[1] < 5*h/6.0 + tol
                return (xin and (yb or yt)) or (yval and xval)

class inside_csc(SubDomain):
        def inside(self,x,on_boundary):
                return (abs(x[0]) < s + tol and (h/6.0 - tol < x[1] < 5*h/6.0 + tol))

class cord(SubDomain):
        def inside(self,x,on_boundary):
                csc = (abs(x[0]) < s + tol and (h/6.0 - tol < x[1] < 5*h/6.0 + tol))
		xval = s-tol < abs(x[0]) < s+tol
                yval = h/6.0-tol < abs(x[1]) < h*5/6.0 + tol
                xin = -s - tol < x[0] < s + tol
                yb = h/6.0 - tol < x[1] < h/6.0 + tol
                yt = 5*h/6.0 - tol < x[1] < 5*h/6.0 + tol
                csc_bnd = (xin and (yb or yt)) or (yval and xval)
                return abs(x[0]) < x0 + tol# and (not csc or csc_bnd)   

class CSF_space(SubDomain):
        def inside(self,x,on_bounary):
                return x0-tol < abs(x[0]) < x1+tol

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


bot_out = bottom_outer()
bot_in = bottom_inner()
out_wall = outer_wall()
in_wall = inner_wall()
top_out = top_outer()
top_in = top_inner()
csc_bnd = csc_boundary()
in_csc = inside_csc()
cord = cord()
CSF_s = CSF_space()

