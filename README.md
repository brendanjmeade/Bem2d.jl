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

![ex_thrusttopo](/docs/src/assets/ex_thrusttopo.png)

## Running inside Docker

Bem2d.jl has its own docker image! Before, continuing, follow these directions to make sure you have the latest version of docker and docker-compose:

* https://docs.docker.com/install/ (also, on Linux, make sure to follow the post-install instructions: https://docs.docker.com/install/linux/linux-postinstall/)
* https://docs.docker.com/compose/install/

To test that out, try running `docker run hello-world`. 

If that works, let's launch a Jupyter notebook with Bem2d.jl installed. Run: `./launch_jupyter` That should first download or build the docker image which includes installing all the dependencies for Bem2d.jl. Afterwards, it'll launch a Jupyter Lab server instance. Find the line in the output saying `[I 19:20:33.891 LabApp] The Jupyter Notebook is running at:`. The next line should tell you what port the server is running on. Open your browser to `localhost:THATPORT` and play around! 
