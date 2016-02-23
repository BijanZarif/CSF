from dolfin import *
import pylab as plt
from pylab import zeros, where, linspace, ones, show, array

C1_Brucker = 10*array([2.26, 1.53, 0.18, 0.19, 0.2, 0.05, -0.39, -0.69, -0.93, -0.66, -0.63, -0.7, -0.23, 2.51])      # 12:00

C1_Brucker2 = 10*array([-0.33, -0.35, -0.51, 0.99, 1.27, 0.83, 0.71, 0.67, -0.15, -0.71, -0.05, -0.21, -0.43, -0.62]) # 6:00


t_Brucker = 1e-3*array([10, 73, 136, 199, 262, 325, 388, 451, 514, 577, 640, 703, 766, 829])
C = C1_Brucker[5:]
C2 = C1_Brucker[:5]
C1_Brucker = plt.append(C,C2)

Cop = C1_Brucker2[5:]
Cop2 = C1_Brucker2[:5]
C1_Brucker2 = plt.append(Cop,Cop2)

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
		idx = plt.find(t_Brucker >= t)[0] - 1

		values[1] = -(yval[idx] + (yval[idx+1]-yval[idx])*(t-tval[idx])/(tval[idx+1]-tval[idx]))
		values[0] = 0

	def value_shape(self):
		return (2,)

