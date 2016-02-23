from dolfin import *
H = 0.41;
L = 2.5;

l = 0.35;
h = 0.02;

xb = 0.6;
yb = 0.19;
Cx = 0.2;
Cy = 0.2;
r = 0.05;
xL = Cy+0.04898979485566357;

eps = 1e-6

class Outlet(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[0],L) and on_bnd

class Inlet(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[0],0.0) and on_bnd

class Top(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],H) and on_bnd

class Bottom(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],0.0) and on_bnd

class Interface(SubDomain):
	def inside(self,x,on_bnd):
		return (near(x[1],yb+h) or near(x[1],yb) or near(x[0],Cx+r+l)) and abs(x[0] - (Cx + r + l/2.0)) < l and abs(x[1] - (yb + 0.5*h)) < h


class Circle(SubDomain):
        def inside(self,x,on_bnd):
                return ((x[0]-Cx)**2 + (x[1]-Cy)**2 <= r**2 + eps) and on_bnd

class Elastic(SubDomain):
	def inside(self,x,on_bnd):
		return abs(x[1]-(yb+h/2)) < h/2 + eps and x[0] > xL - eps and x[0] < xb + eps

class Elastic_to_Circle(SubDomain):
	def inside(self,x,on_bnd):
		E = Elastic()
		C = Circle()
		return C.inside(x,on_bnd) and E.inside(x,on_bnd)		


outlet = Outlet()
inlet = Inlet()
top = Top()
bottom = Bottom()
interface = Interface()
circle = Circle()

def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)
