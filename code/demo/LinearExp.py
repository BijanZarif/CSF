from dolfin import *
import pylab as plt

class MyExpression0(Expression):
	def __init__(self,t_Brucker,C1_Brucker,t):
		self.t = t
		self.t_Brucker = t_Brucker
		self.C1_Brucker = C1_Brucker

	def eval(self,values,x):
		print 'evaluation', values, x, self.t
		t = self.t
		t_Brucker = self.t_Brucker
		C1_Brucker = self.C1_Brucker
		while float(t) > self.t_Brucker[-1]:
			t -= self.t_Brucker[-1]
		tval = t_Brucker
		yval = C1_Brucker
		idx = plt.find(t_Brucker >= float(t))[0] - 1
		values[0] = -(yval[idx] + (yval[idx+1]-yval[idx])*(t-tval[idx])/(tval[idx+1]-tval[idx]))
