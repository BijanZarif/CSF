#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *
from cbcpost import PostProcessor

from os import path
#set_log_active(False)

files = [path.join(path.dirname(path.realpath(__file__)),"cbcflow-data/cylinder_0.6k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"cbcflow-data/cylinder_2k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"cbcflow-data/cylinder_8k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"cbcflow-data/cylinder_32k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"cbcflow-data/cylinder_129k.xml.gz"),
        ]

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 10.0)

class Wall(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and not (near(x[0],10.0) or near(x[0], 0.0))


class FlowAroundCylinder(NSProblem):
    "Flow around a cylinder in 2D."

    @classmethod
    def default_params(cls):
        "Default parameters overwriting and adding to NSProblem.default_params()"
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=3.0,
            dt=0.005,
            # Physical parameters
            rho=1.0,
            mu=1./100#2.414*10**-5*10**(247.8/(310-140)),
            )
        params.update(
            # Spatial parameters
            refinement_level=0,
            )
        return params

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Load mesh
        mesh = Mesh('spinal_cord.xml')#Mesh(files[self.params.refinement_level])
        mesh = refine(mesh)
        # The mesh can also be generated on the fly. This has been
        # commented out because the mesh generator is non-deterministic and thus
        # unsuitable for the test suites.
        """
        refinement_levels=[32,64,128,256,512]
        N = refinement_levels[self.params.refinement_level]
        # Create mesh
        r = Rectangle(0,0, 10, 1)
        c = Circle(2.0, 0.5, 0.12)
        mesh = Mesh(r-c, N)
        """

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(3)
        Wall().mark(facet_domains, 0)
        #Cylinder().mark(facet_domains, 0)
        LeftBoundary().mark(facet_domains, 1)
        RightBoundary().mark(facet_domains, 2)
        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)
        

    def initial_conditions(self, spaces, controls):
        "Setting the flow at rest as initial conditions"
        c0 = Constant(0)
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):

        c0 = Constant(0)
        c1 = Constant(1)
        p0 = Constant(0)
	p1 = Expression('-10.0*cos(2*t)',t=t)
	p2 = Expression('100.0*sin(2*pi*t)',t=t)

        # Create inflow and no-slip boundary conditions for velocity
        u0 = Expression(('2*cos(4*pi*t)','0.0'),t=t)
        
        inflow = (u0, 1)
        noslip = ([c0, c0], 0)

        # Create boundary conditions for pressure
        bcpR = (c0, 2)
        bcpR2 = (p2, 2)
	bcpL = (p1, 1)
	
        # Collect and return
        bcu = [noslip]
        bcp = [bcpR2,bcpL]
        return (bcu, bcp)


def main():
    # Create problem and scheme instances
    problem = FlowAroundCylinder({"refinement_level": 2})
    scheme = IPCS_Stable()

    # Create postprocessor instance pointing to a case directory
    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    postprocessor = PostProcessor({"casedir": casedir})

    # Creating fields to plot and save
    plot_and_save = dict(plot=False, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        ]

    # Add fields to postprocessor
    postprocessor.add_fields(fields)

    # Create NSSolver instance and solve problem
    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()


if __name__ == "__main__":
    main()
