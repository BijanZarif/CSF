from dolfin import *
from mshr import *

h = 60.0

x0 = 5.0
x1 = 9.0
s = 3.0

P0 = Point(0,0,0)
P1 = Point(0,0,h/6.0)
P2 = Point(0,0,5*h/6.0)
P3 = Point(0,0,h)
P4 = Point(0,0,h/2)

cyst = Ellipsoid(P4,s,s,2*h/6)

solid = Cylinder(P0,P3,x0,x0,40)-cyst

fluid = Cylinder(P0,P3,x1,x1,40)-solid-cyst

res = 20.0



a = CSGCGALDomain3D(cyst)
b = CSGCGALDomain3D(solid)
c = CSGCGALDomain3D(fluid)


generator = TetgenMeshGenerator3D()
generator.parameters["preserve_surface"] = True
generator.parameters["mesh_resolution"] = res

mesh_c = generator.generate(a)
mesh_s = generator.generate(b)
mesh_f = generator.generate(c)
'''
inner = DolfinMeshUtils.merge_meshes(mesh_s, mesh_c)
outer = DolfinMeshUtils.merge_meshes(mesh_s, mesh_f)
mesh = DolfinMeshUtils.merge_meshes(inner, outer)


mesh_c = generate_mesh(cyst,res)
mesh_s = generate_mesh(solid,res)
mesh_f = generate_mesh(fluid,res)
'''
f0 = File('RESULTS/cyst.xml')
f1 = File('RESULTS/solid.xml')
f2 = File('RESULTS/fluid.xml')
f0 << mesh_c
f1 << mesh_s
f2 << mesh_f

#F = File('RESULTS/merged.xml.gz')
#F << mesh



