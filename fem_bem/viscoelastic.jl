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
VE_z = -16e3 # Depth below which we use a Maxwell rheology
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
xgrid, ygrid = obsgrid(L + offset, B + offset, R - offset, T - offset, npts);

### Slip from fault only
T_TB_F, H_TB_F = PUTC(slip2dispstress, elsbox, BRTL_idx, idx["F"], mu, nu)
Fslip = [100; 100]; # y-direction slip only
Uslip = T_TB_F * Fslip
Tslip = H_TB_F * Fslip;

# Kernels and solve
T_B_BRTL, H_B_BRTL = PUTC(slip2dispstress, elsbox, idx["B"], BRTL_idx, mu, nu)
T_RTL_BRTL, H_RTL_BRTL = PUTC(slip2dispstress, elsbox, RTL_idx, BRTL_idx, mu, nu)
bcs_fault_only = zeros(8 * nside)
bcs_fault_only[1:2 * nside] = -Uslip[1:2 * nside] # Bottom
bcs_fault_only[2 * nside + 1:4 * nside] = -Tslip[2 * nside + 1:4 * nside] # Right
bcs_fault_only[4 * nside + 1:6 * nside] = -Tslip[4 * nside + 1:6 * nside] # Top
bcs_fault_only[6 * nside + 1:8 * nside] = -Tslip[6 * nside + 1:8 * nside] # Left
THbox = [T_B_BRTL ; H_RTL_BRTL];

THboxinv = inv(THbox);

Ueff_fault_only = THbox \ bcs_fault_only
UTB, STB = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, BRTL_idx,
                           Ueff_fault_only[1:2:end], Ueff_fault_only[2:2:end], mu, nu)
UF, SF = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, idx["F"],
                         Fslip[1:2:end], Fslip[2:2:end], mu, nu)
Ufaultonly = UTB .+ UF
Sfaultonly = STB .+ SF
twopanel(xgrid, ygrid, npts, Ufaultonly, Sfaultonly, idx, elsbox, figsize=(12, 6), ylim=[B,T]);

fenics_mesh = RectangleMesh(Point((L - 5000, T + 5000)), Point((R + 5000, B - 5000)), 60, 60)
scalar_fs = FunctionSpace(fenics_mesh, "P", 1)
vector_fs = VectorFunctionSpace(fenics_mesh, "P", 2)
tensor_fs = FunctionSpace(FEniCS.fenics.TensorFunctionSpace(fenics_mesh.pyobject, "P", 1))

include("Helpers.jl")

x_dofs = py"lambda x: x.tabulate_dof_coordinates()"(scalar_fs.pyobject);
_, stress_from_fault = constdispstress(slip2dispstress, x_dofs[:,1], x_dofs[:,2], elsbox, idx["F"],
                                       Fslip[1:2:end], Fslip[2:2:end], mu, nu)

# no viscous behavior in the surface layer, mu/eta below VE_z
# We can't use a threshold for eta because we divide by eta and dividing by zero is bad.
fenics_mu_over_eta = interpolate(Expression("0.0 + mu/eta*(x[1] <= VE_z)", degree=0, mu=mu, eta=eta, VE_z=VE_z), scalar_fs)
fenics_Ve = FeFunction(vector_fs)
fenics_disp = FeFunction(vector_fs)

v = TestFunction(vector_fs)
bc = DirichletBC(vector_fs, Constant((0, 0)), "on_boundary")
function fenics_lhs_assemble(V)
    u = TrialFunction(V)
    a = inner(fenics_stress(u, lambda, mu), fenics_strain(v)) * dx
    A = assemble(a)
    apply(bc, A)
    A
end
function fenics_rhs_assemble(V, f)
    L = FEniCS.dot(f, v) * dx
    b = assemble(L)
    apply(bc, b)
    b
end
lhs = fenics_lhs_assemble(vector_fs)
lu_solver = py"""
from fenics import *
def create_lu_solver(lhs):
    s = LUSolver(lhs, method='umfpack')
    s.parameters['symmetric'] = True
    return s
"""
lu_solver = py"create_lu_solver"(lhs.pyobject)
function fenics_lu_solve(lu_solver, u, rhs)
    py"lambda lu_solver, u, rhs: lu_solver.solve(u.vector(), rhs)"(lu_solver, u.pyobject, rhs.pyobject)
end

# Step 1: solve FEM+BEM problem
function solve_elastic_problem(dofs)
    vector(fenics_Ve).set_local(dofs)
    rhs = fenics_rhs_assemble(vector_fs, fenics_Ve)
    fenics_lu_solve(lu_solver, fenics_disp, rhs)
        
    fenics_sigma = fenics_stress(fenics_disp, lambda, mu)
    sxx_sxy, sxy_syy = fenics_tensor_components(vector_fs, fenics_sigma)

    bcs_B = -fenics_eval_u(fenics_disp, idx["B"], elsbox.xcenter, elsbox.ycenter)
    bcs_RTL_S = -fenics_eval_s(sxx_sxy, sxy_syy, RTL_idx, elsbox.xcenter, elsbox.ycenter)
    bcs_RTL = stress_to_trac(bcs_RTL_S, RTL_idx, elsbox.xnormal, elsbox.ynormal);

    bcs_ve_fem = zeros(2 * elsbox.endidx - 2)
    bcs_ve_fem[1:2 * nside] = bcs_B
    bcs_ve_fem[2 * nside + 1:end] = bcs_RTL;

    bcscombined = bcs_fault_only .+ bcs_ve_fem    
    Ueffcombined = THbox \ bcscombined
    
    fenics_sigma, Ueffcombined
end

# Step 2: convert BEM stress to a fenics object
# Three sources of stress: u_eff, fault, FEM
function calc_total_stress(fenics_sigma, Ueffcombined)
    _, stress_from_ueff = constdispstress(slip2dispstress, x_dofs[:,1], x_dofs[:,2], elsbox, BRTL_idx,
                                          Ueffcombined[1:2:end], Ueffcombined[2:2:end], mu, nu)
    total_bem_stress = stress_from_ueff .+ stress_from_fault
    
    # set stress values to zero when they are outside the box domain. TODO: I don't think this is right solution even though it's easy.
    total_bem_stress[x_dofs[:,1] .<= L, :] .= 0
    total_bem_stress[x_dofs[:,1] .>= R, :] .= 0
    total_bem_stress[x_dofs[:,2] .<= B, :] .= 0
    total_bem_stress[x_dofs[:,2] .>= T, :] .= 0

    fenics_stress_dofs = zeros(4, size(x_dofs)[1])
    fenics_stress_dofs[1,:] = total_bem_stress[:,1]
    fenics_stress_dofs[2,:] = total_bem_stress[:,3]
    fenics_stress_dofs[3,:] = total_bem_stress[:,3]
    fenics_stress_dofs[4,:] = total_bem_stress[:,2];
    bem_sigma = FeFunction(tensor_fs)
    vector(bem_sigma).set_local(reshape(fenics_stress_dofs, (:)))
    total_sigma = bem_sigma + fenics_sigma
end

# Step 3: calculate the stress deviator and dVdt. Extract dofs from fenics function and return the array.
function stress_to_dVdt(total_sigma)
    # https://en.wikipedia.org/wiki/Cauchy_stress_tensor#Stress_deviator_tensor
    stress_deviator = total_sigma - (1.0 / 3) * FEniCS.tr(total_sigma) * Identity(2);
    dVdt = project(fenics_mu_over_eta * div(stress_deviator), vector_fs)
    vector(dVdt).get_local()
end

function volume_disp_stress(Ueffcombined)
    UTB, STB = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, BRTL_idx,
                               Ueffcombined[1:2:end], Ueffcombined[2:2:end], mu, nu)
            
    fenics_sigma = fenics_stress(fenics_disp, lambda, mu)
    sxx_sxy, sxy_syy = fenics_tensor_components(vector_fs, fenics_sigma)
    
    Ufem = transpose(reshape(fenics_eval_u(fenics_disp, 1:length(xgrid), xgrid, ygrid), (2, :)))
    Sfem = transpose(reshape(fenics_eval_s(sxx_sxy, sxy_syy, 1:length(xgrid), xgrid, ygrid), (3, :)));
    
    Utotal = UTB .+ UF .+ Ufem
    Stotal = STB .+ SF .+ Sfem
    Utotal, Stotal
end

function calc_dVdt(Ve_force_dofs)
    fenics_sigma, Ueffcombined = solve_elastic_problem(Ve_force_dofs)
    total_sigma = calc_total_stress(fenics_sigma, Ueffcombined)
    dVdt  = stress_to_dVdt(total_sigma)
    dVdt, Ueffcombined
end

function main()
    V = zeros(vector(fenics_Ve).size())
    dt = siay * 2
    t = 0
    Uhistory = []
    for i in 1:100
        println(t)
        dVdt, Ueffcombined = calc_dVdt(V)
        V += dt .* dVdt
        t += dt
        Utotal, Stotal = volume_disp_stress(Ueffcombined)
        append!(Uhistory, [Utotal])
    # @time twopanel(xgrid, ygrid, npts, Utotal, Stotal, idx, elsbox, figsize=(12,6), ylim=[B,T]);
    end
    Uhistory = reshape(hcat(Uhistory...), (:, 2, length(Uhistory)));

    # sanity check: velocity measured in mm/yr
    Vhistory = (Uhistory[:,:,2:end] - Uhistory[:,:,1:end - 1]) ./ dt * siay * 1000;

    maxV = dropdims(maximum(sqrt.(Vhistory[:,1,:].^2 + Vhistory[:,2,:].^2), dims=(1)), dims=(1));

    maxV

    plot(maxV)

    for i in 1:10:41
        figure()
        field = sqrt.(Vhistory[:, 1, i].^2 + Vhistory[:, 2, i].^2)
        ncontours = 10
        lowfield = minimum(field)
        highfield = maximum(field)
        ncontours = LinRange(lowfield, highfield, 11)    

        xlim = [-20000 20000]
        scale = 1.0
        fieldmax = maximum(@.abs(field))
        contourf(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
             reshape(field, npts, npts), ncontours,
             vmin=lowfield, vmax=highfield,
             cmap=get_cmap("magma"))
        clim(lowfield, highfield)
        colorbar(fraction=0.020, pad=0.05, extend="both", label=L"$||u||$ (m)")
        contour(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
            reshape(field, npts, npts), ncontours,
            vmin=lowfield, vmax=highfield,
            linewidths=0.25, colors="w")
        gca().set_aspect("equal")
        gca().set_xlim([-20000, 20000])
        gca().set_ylim([-50000, 0])
        xlabel(L"$x$ (m)")
        ylabel(L"$y$ (m)")
    end

    plot(log10.(maxV))

end
main()


