from matplotlib.pylab import *




a = 8.0
b = 0.5
r = sqrt(a**2+b**2)
x = linspace(-r*sqrt(b),r*sqrt(b),101)
y1 = sqrt(r**2-x**2/b)*sqrt(a)
y2 = -y1


plot(x,y1,x,y2)
axis([-20, 20, -20, 20])
show()
