from dolfin import *
import os

os.system('gmsh -3 drosdal3D.geo')
os.system('dolfin-convert drosdal3D.msh drosdal3D.xml')



mesh = Mesh("drosdal3D.xml")                                           
subdomains = MeshFunction("size_t", mesh, "drosdal3D_physical_region.xml")    
boundaries = MeshFunction("size_t", mesh, "drosdal3D_facet_region.xml")
hdf = HDF5File(mesh.mpi_comm(), "drosdal3D.h5", "w")
hdf.write(mesh, "/mesh")
hdf.write(subdomains, "/subdomains")
hdf.write(boundaries, "/boundaries")
