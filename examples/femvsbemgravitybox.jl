using FEniCS
using PyPlot


# Calculate strain from displacment
function epsilon(u)
     0.5 * (nabla_grad(u) + Transpose(nabla_grad(u)))
end

# Calculate stress from displacement
# s = sigma(u) - (1.0/3.0) * tr(sigma(u)) * Identity(d)
function sigma(u)
     lambda * nabla_div(u) * Identity(d) + 2 * mu * epsilon(u)
end


function femvsbemgravitybox()
    mu = 1
    lambda = 1
    rho = 1
    g = 9.81

    # Create mesh and define function space
    mesh = RectangleMesh(Point((0.0, 0.0)), Point((1.0, 1.0)), 10, 10)
    V = VectorFunctionSpace(mesh, "P", 1)
    c = Constant((0, 0))
    bc = DirichletBC(V, c, "on_boundary && x[1]<1E-14") # What BCs are bing set???

    # Solve
    u = TrialFunction(V)
    d = geometric_dimension(u) # space dimension
    v = TestFunction(V)
    f = Constant((0, -rho*g)) # Vector of uniform body force
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
