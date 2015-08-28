from dolfin import *
import numpy as np
set_log_active(False)
mesh = Mesh('cord_w_csc.xml')

#mesh = refine(mesh)
#mesh = refine(mesh)	

rho = 1.0
nu = 10**-3
h = 0.06
x0 = 0.005
x1 = 0.009
s = 0.001

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
                yval = 0.01-tol < abs(x[1]) < 0.05 + tol
                xin = -s - tol < x[0] < s + tol
                yb = 0.01 - tol < x[1] < 0.01 + tol
                yt = 0.05 - tol < x[1] < 0.05 + tol
                return (xin and (yb or yt)) or (yval and xval)

class inside_csc(SubDomain):
        def inside(self,x,on_boundary):
                return (abs(x[0]) < s + tol and (0.01 - tol < x[1] < 0.05 + tol))

class cord(SubDomain):
        def inside(self,x,on_boundary):
                csc = (abs(x[0]) < s + tol and (0.01 - tol < x[1] < 0.05 + tol))
                return abs(x[0]) < x0 and not csc

class CSF_space(SubDomain):
        def inside(self,x,on_bounary):
                return 0.005-tol < abs(x[0]) < 0.009+tol
