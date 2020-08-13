# Bem2d.jl
Bem2d.jl is a Julia library for two-dimensional plane strain linear elastic boundary element (BEM) calculations related to earthquakes and tectonics.

## Package Features
  - Analytic quadratic slip and traction elements (Portella et al., 1992)
  - Surface topography & non-planar fault surfaces
  - Unidirectional gravity (Pape and Banerjee, 1987)
  - Kinematically consistent quasi-dynamic earthquake cycles with slip and aging laws

## Basic BEM calculation
A simple starting example of a thrust fault with a topographic surface can be generated with:

```Julia
using Bem2d
include(“examples/thrusttopo.jl”)
```
should produce this figure:

![thrusttopo](/docs/src/assets/ex_thrusttopo.png)

## Examples
These serve as reference examples (and some as effectively tests) with the intent of making clear how to build simple to complex BEM models

`bemvsdeepokada.jl` - [works on Linux, because of okada_wrapper dependency] - 
Comparision between fundamental fullspace 2D elastic disloation soluitons and Okada's half space 3D solutions. 

`braziltest.jl` - [working] - Traction only compression of a homogeneous disc (Brazil test)

`constvsquadelements.jl` - [working] - Basic comparision of constant slip vs. 3 node quadratic elements.

`discmaterial.jl` - [working] - A simple two domain problem.  An inner annulus surrounded by a an unbounded region with different material properties.  Crouch and Starfield solve this for using fictious force vs. the DDM approach here.

`dislocationinabox.jl` - [not working] - Single fault inside of a box with traction free boundary conditions on the sides and top with a no slip boundary at the bottom.

`dislocationinaboxgravity.jl` - [not working] - Single fault inside of a box with gravity. Traction free boundary conditions on the sides and top with a no slip boundary at the bottom.

`gravitydislocation.jl` - [not working] - A box with gravity and fault with a known slip.

`gravitysquare.jl` - [working] - Homogenous square compressing under the force of gravity with a no slip bottom boundary condition.

`kelvin.jl` - [working] - Comparison of point source kelvin solution with finite source DDM approximation.

`layeredboxfault.jl` - [not working] - Box with a fault and two domains with differing material properties.

`nearsurfaceerror.jl` - [works on Linux, because of okada_wrapper dependency] - Comparision of constant element vs. 3 node quadratic element errors when compared with Okada.

`okadaindirect.jl` - [works on Linux, because of okada_wrapper dependency] - Comparision of Okada's half space solution with a BEM implementation.

`qdplanar.jl` - [not working] - Planar quasi-dynamic rupture example.

`qdtopography.jl` - [not working] - Quasi-dynamic rupture example with topography.

`thrustfaultfreesurface.jl` - [working] - Older version of Okada comparision at free surface.  Needs to updated.

`thrustfaulttopography.jl` - [working] - A thrust fault with a traction free topographic surface

`weldedcircle.jl` - [working] - A basic two domain example consisting of a circle cut into upper and lower halves with differing material properties.  Traction boundary conditions apply on the top of the circle and no-slip boundary conditions on the bottom half.

`weldedcirclefault.jl` - [working] - A basic two domain example consisting of a circle cut into upper and lower halves with differing material properties and fault in the upper domain. The upper half of the circle boundary is traction and no-slip boundary conditions on the bottom half.


## Running inside Docker

Bem2d.jl has its own docker image! Before, continuing, follow these directions to make sure you have the latest version of docker and docker-compose:

* https://docs.docker.com/install/ (also, on Linux, make sure to follow the post-install instructions: https://docs.docker.com/install/linux/linux-postinstall/)
* https://docs.docker.com/compose/install/

To test that out, try running `docker run hello-world`. 

If that works, let's launch a Jupyter notebook with Bem2d.jl installed. Run: `./launch_jupyter` That should first download or build the docker image which includes installing all the dependencies for Bem2d.jl. Afterwards, it'll launch a Jupyter Lab server instance. Find the line in the output saying `[I 19:20:33.891 LabApp] The Jupyter Notebook is running at:`. The next line should tell you what port the server is running on. Open your browser to `localhost:THATPORT` and play around! 
