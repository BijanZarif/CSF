from matplotlib.pylab import *

x = linspace(0,1,6)

p0 = [1,0]
p = [0,1,0]

plot(x[0:2],p0)
for i in range(2):
	plot(x[i:i+3],p)
l = ['$\phi_0$','$\phi_1$','$\phi_2$']
legend(l)
plot(x,x*0,'o')
plot(x[:2],2*[0],'r')
plot(x[1:3],2*[0],'b')
plot(x[2:4],2*[0],'g')
plot(x[3:],3*[0],'r')
axis([0, 1, -0.1, 1.1])
show()
