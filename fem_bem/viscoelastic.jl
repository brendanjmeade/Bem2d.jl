using Revise
using PyPlot
using Infiltrator
using LinearAlgebra
using Bem2d
using PyCall
using FEniCS

include("Helpers.jl")

mu = 30e9
lambda = 30e9
eta = 1e19 # just made up a number!
nu = 0.25
g = 9.81
rho = 2700
siay = 3600 * 24 * 365.2425

# Element geometries and data structures 
elsbox = Elements(Int(1e5))

# Box geometry
B = -50e3 # Bottom
R = 20e3 # Right
T = 0e3 # Top
L = -20e3; # Left
VE_z = -20e3 # Depth below which we use a Maxwell rheology
nside = 50 # number of elements per box side

addelsez!(elsbox, discretizedline(L, B, R, B, nside)..., "B")
addelsez!(elsbox, discretizedline(R, B, R, T, nside)..., "R")
addelsez!(elsbox, discretizedline(R, T, L, T, nside)..., "T")
addelsez!(elsbox, discretizedline(L, T, L, B, nside)..., "L")

# Fault geometry
nfault = 1
x1, y1, x2, y2 = discretizedline(-10e3, -10e3, -5e3, -5e3, nfault) # 45 degree dipping fault
addelsez!(elsbox, x1, y1, x2, y2, "F")

idx = getidxdict(elsbox);
BRTL_idx = [idx["B"] ; idx["R"] ; idx["T"] ; idx["L"]]
RTL_idx = [idx["R"] ; idx["T"] ; idx["L"]];

npts = 50
offset = 1
xgrid, ygrid = obsgrid(L+offset, B+offset, R-offset, T-offset, npts);

### Slip from fault only
T_TB_F, H_TB_F = PUTC(slip2dispstress, elsbox, BRTL_idx, idx["F"], mu, nu)
Fslip = [100; 100]; # y-direction slip only
Uslip = T_TB_F * Fslip
Tslip = H_TB_F * Fslip;

# Kernels and solve
T_B_BRTL, H_B_BRTL = PUTC(slip2dispstress, elsbox, idx["B"], BRTL_idx, mu, nu)
T_RTL_BRTL, H_RTL_BRTL = PUTC(slip2dispstress, elsbox, RTL_idx, BRTL_idx, mu, nu)
bcs_fault_only = zeros(8 * nside)
bcs_fault_only[1:2*nside] = -Uslip[1:2*nside] # Bottom
bcs_fault_only[2*nside+1:4*nside] = -Tslip[2*nside+1:4*nside] # Right
bcs_fault_only[4*nside+1:6*nside] = -Tslip[4*nside+1:6*nside] # Top
bcs_fault_only[6*nside+1:8*nside] = -Tslip[6*nside+1:8*nside] # Left
THbox = [T_B_BRTL ; H_RTL_BRTL];

THboxinv = inv(THbox);

Ueff_fault_only = THbox \ bcs_fault_only
UTB, STB = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, BRTL_idx,
                           Ueff_fault_only[1:2:end], Ueff_fault_only[2:2:end], mu, nu)
UF, SF = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, idx["F"],
                         Fslip[1:2:end], Fslip[2:2:end], mu, nu)
Ufaultonly = UTB .+ UF
Sfaultonly = STB .+ SF
twopanel(xgrid, ygrid, npts, Ufaultonly, Sfaultonly, idx, elsbox, figsize=(12,6), ylim=[B,T]);

fenics_mesh = RectangleMesh(Point((L - 5000, T + 5000)), Point((R + 5000, B - 5000)), 100, 100)
scalar_fs = FunctionSpace(fenics_mesh, "P", 1)
vector_fs = VectorFunctionSpace(fenics_mesh, "P", 2) # TODO: can i make this order = 1?
tensor_fs = FunctionSpace(FEniCS.fenics.TensorFunctionSpace(fenics_mesh.pyobject, "P", 1))

Ueff = Ueff_fault_only;
x_dofs = py"lambda x: x.tabulate_dof_coordinates()"(scalar_fs.pyobject);

_, stress_from_ueff = constdispstress(slip2dispstress, x_dofs[:,1], x_dofs[:,2], elsbox, BRTL_idx,
                                      Ueff[1:2:end], Ueff[2:2:end], mu, nu)

# set stress_from_ueff values to zero when they are outside the box domain.
stress_from_ueff[x_dofs[:,1] .<= L, :] .= 0
stress_from_ueff[x_dofs[:,1] .>= R, :] .= 0
stress_from_ueff[x_dofs[:,2] .<= B, :] .= 0
stress_from_ueff[x_dofs[:,2] .>= T, :] .= 0

fenics_stress_dofs = zeros(4,size(x_dofs)[1])
fenics_stress_dofs[1,:] = stress_from_ueff[:,1]
fenics_stress_dofs[2,:] = stress_from_ueff[:,3]
fenics_stress_dofs[3,:] = stress_from_ueff[:,3]
fenics_stress_dofs[4,:] = stress_from_ueff[:,2];
sigma = FeFunction(tensor_fs)
py"lambda x, vs: x.vector().set_local(vs)"(sigma.pyobject, reshape(fenics_stress_dofs, (:)))

fenics_plot_vec2d(FEniCS.dot(sigma, Constant((1.0, 0.0))))
fenics_plot_vec2d(FEniCS.dot(sigma, Constant((0.0, 1.0))))

# no viscous behavior in the surface layer, mu/eta below VE_z
# We can't use a threshold for eta because we divide by eta and dividing by zero is bad.
fenics_mu_over_eta = interpolate(Expression("0.0 + mu/eta*(x[1] <= VE_z)", degree=0, mu=mu, eta=eta, VE_z=VE_z), scalar_fs)

# https://en.wikipedia.org/wiki/Cauchy_stress_tensor#Stress_deviator_tensor
stress_deviator = sigma - (1.0/3)*FEniCS.tr(sigma)*Identity(2);

dVdt = project(fenics_mu_over_eta * div(stress_deviator), vector_fs)

fenics_plot_vec2d(dVdt)

function calc_dVdt()
    # Step 1: calculate BEM stress
    # Step 2: convert BEM stress to a fenics object
    # Step 3: calculate the stress deviator
    # Step 4: calculate dVdt from the stress deviator.
end
