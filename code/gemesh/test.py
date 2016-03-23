from dolfin import *



c = Cone(Point(0, 0, -1), Point(0, 0, 1), .5, 1.0)
m = Mesh(c, 16)

plot(m)
interactive()
