from dolfin import *

mesh = UnitSquareMesh(10)

V = FunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)

VQ = V*Q

up = TrialFunction(VQ)
v,q = TestFunctions(VQ)

noslip = Constant(('0.0'),('0.0'))

bcu = DirichletBC(VQ.sub(0),
