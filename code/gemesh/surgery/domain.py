from dolfin import *

L = 30;
H = 6;
r = 4;
R = 4;

s = 10;

eps = 1e-8

class Solid(SubDomain):
	def inside(self,x,on_bnd):
		return x[1] > H - eps

class SolidNoSlip(SubDomain):
	def inside(self,x,on_bnd):
		return x[1] > H - eps and on_bnd

class RigidInterface(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],H) and abs(x[0]-(s+r)) > r - eps

class SoftInterface(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[1],H) and abs(x[0]-(s+r)) < r + eps

class InletF(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[0],0) and x[1] < H + eps

class Bottom(SubDomain):
	def inside(self,x,on_bnd):
		return (near(x[1],0) or (x[1]/R)**2 + ((x[0]-(s+r))/r)**2 < 1+eps) and on_bnd
class Outlet(SubDomain):
	def inside(self,x,on_bnd):
		return near(x[0],L) and on_bnd
