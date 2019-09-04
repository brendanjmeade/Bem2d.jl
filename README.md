# Bem2d.jl
Bem2d.jl is a Julia library for two-dimensional boundary element (BEM) calculations related to earthquakes and tectonics.

## Package Features
  - Constant and quadratic slip and traction elastic elements.
  - Surface topography & non-planar fault surfaces
  - Quasi-dynamic earthquake cycles with both slip and aging laws.

## Basic BEM calculation
A simple starting example of a thrust fault with a topographic surface can be generated with:

```Julia
using Bem2d
include(“ex_thrusttopo.jl”)
```
should produce this figure:

![ex_thrusttopo](/docs/src/assets/ex_thrusttopo.png)
