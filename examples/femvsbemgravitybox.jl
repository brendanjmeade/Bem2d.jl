using FEniCS
using PyPlot

function femvsbemgravitybox()
    mu = 3e10
    lambda = mu
    rho = 2700
    g = 9.81
    L = 10e3

    # Create mesh and define function space
    mesh = RectangleMesh(Point((-L, 0)), Point((L, 2 * L)), 100, 100)
    V = VectorFunctionSpace(mesh, "P", 1)
    c = Constant((0, 0))
    bc = DirichletBC(V, c, "on_boundary && x[1]<1E-14") # What BCs are bing set???

    # Solve
    u = TrialFunction(V)
    d = geometric_dimension(u) # space dimension
    v = TestFunction(V)
    f = Constant((0, -rho * g)) # Vector of uniform body force
    T = Constant((0, 0))
    a = inner(sigma(u), epsilon(v)) * dx
    L = dot(f, v) * dx + dot(T, v) * ds
    u = FeFunction(V)
    lvsolve(a, L, u, bc)

    # Extract plottable quantities
    V = FunctionSpace(mesh, "P", 1)
    u_magnitude = sqrt(dot(u, u))
    u_magnitude = project(u_magnitude, V)
    close("all")
    plothandle = FEniCS.Plot(u_magnitude)
    colorbar(plothandle)
end
femvsbemgravitybox()
