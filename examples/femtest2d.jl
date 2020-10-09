using FEniCS
using PyPlot

mu = 1
rho = 1
# gamma = 0.4*delta^2
# beta = 1.25
lambda = 1
# g = gamma
g = 9.81


# Create mesh and define function space
mesh = UnitSquareMesh(20, 20)
V = VectorFunctionSpace(mesh, "P", 1)
c = Constant((0, 0))
bc = DirichletBC(V, c, "on_boundary && x[1]<1E-14")
tol = 1E-14

# Define strain
function epsilon(u)
     0.5*(nabla_grad(u)+Transpose(nabla_grad(u)))
end

# Define stress
function sigma(u)
     lambda * nabla_div(u) * Identity(d) + 2 * mu * epsilon(u)
end

u = TrialFunction(V)
d = geometric_dimension(u) # space dimension
v = TestFunction(V)
f = Constant((0, -rho*g)) # Vector of uniform body force
T = Constant((0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds
u = FeFunction(V)
lvsolve(a, L, u, bc)

s = sigma(u) - (1.0/3.0) * tr(sigma(u)) * Identity(d)
V = FunctionSpace(mesh, "P", 1)
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)

FEniCS.Plot(u_magnitude)
