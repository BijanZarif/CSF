from sympy import *

x3,C = symbols('x3 C',real=True)
x,h,t,w,nu = symbols('x h t w nu', real=True, positive=True)
k = sqrt(w/(2*nu))

def cc(x):
	return cos(x)*cosh(x)
def ss(x):
	return sin(x)*sinh(x)

f1 = cc(k*x3)*cc(k*h) + ss(k*x3)*ss(k*h)
f2 = cc(k*x3)*ss(k*h) - ss(k*x3)*cc(k*h)
f3 = cc(w)**2 + ss(w)**2 

u = -C/w*((1-f1/f3)*sin(w*t) + f2/f3*cos(w*t))
'''
u1 = re(I*C/w*(1-cosh(sqrt(I*w/nu)*x3)/cosh(sqrt(I*w/nu)*h)))

u = re(u1*exp(I*w*t))

'''

d2u = nu*diff(diff(u,x3),x3)
dt = diff(u,t)

print simplify(d2u-dt)




