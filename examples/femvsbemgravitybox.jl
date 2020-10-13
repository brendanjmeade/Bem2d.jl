using Revise
using FEniCS
using PyPlot
using Infiltrator
using Bem2d

# Calculate strain
function strain(u)
    0.5 * (nabla_grad(u) + Transpose(nabla_grad(u)))
end


# Calculate
function stress(u, lambda, mu)
     lambda * nabla_div(u) * Identity(2) + 2 * mu * strain(u)
end


function femvsbemgravitybox()
    mu = 3e10
    lambda = mu
    rho = 2700
    g = 9.81
    width = 10e3

    ###
    ### FEM solution
    ###
    # Create mesh and define function space
    mesh = RectangleMesh(Point((-width, 0)), Point((width, 2 * width)), 10, 10)
    V = VectorFunctionSpace(mesh, "P", 1)
    bc = DirichletBC(V, Constant((0, 0)), "on_boundary && x[1]<1E-14") # What BCs are bing set???
    
    # Solve
    u = TrialFunction(V)
    d = geometric_dimension(u) # space dimension
    v = TestFunction(V)
    f = Constant((0, -rho * g)) # Vector of uniform body force
    T = Constant((0, 0))
    a = inner(stress(u, lambda, mu), strain(v)) * dx
    L = dot(f, v) * dx + dot(T, v) * ds
    u = FeFunction(V)
    lvsolve(a, L, u, bc)

    # Extract plottable quantities
    V = FunctionSpace(mesh, "P", 1)
    u_magnitude = sqrt(dot(u, u))
    u_magnitude = project(u_magnitude, V)
    close("all")
    figure(figsize=(10, 5))
    subplot(1, 2, 1)
    plothandle = FEniCS.Plot(u_magnitude)
    colorbar(plothandle)
    title("FEM solution")

    # Evaluating the solution at some other point
    # I can extract array values with get_array(u)
    # I should be able to get vertex coords with: map = V.dofmap().vertex_to_dof_map(mesh)
    # ...but that errors out.
    # interppoint = Point((8000, 8000))
    # interpval = u(interppoint)

    figure()
    plot(get_array(u)[2:2:end], "r.")

    @infiltrate

    return
    

    ###
    ### BEM solution
    ### 
    alpha = 7e-8 # scalar preconditioner for traction terms
    nu = 0.25
    nels = 20
    npts = 100
    offset = 1
    x, y = obsgrid(-width+offset, -2*width+offset, width-offset, 0-offset, npts)

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-width, -2*width, width, -2*width, nels) # Bottom
    addelsez!(els, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(width, -2*width, width, 0, nels) # Right hand side
    addelsez!(els, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(width, 0, -width, 0, nels) # Top
    addelsez!(els, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-width, 0, -width, -2*width, nels) # Left hand side
    addelsez!(els, x1, y1, x2, y2, "L")
    
    # Common indexing
    idx = getidxdict(els) # Should this return "all" - YES TODO
    bcidxU = idx["B"] # Boundaries with *displacement* BCs
    bcidxT = [idx["R"] ; idx["T"]; idx["L"]] # Boundaries with *traction* BCs
    bcidxall = collect(1:1:els.endidx) # All boundaries

    # Gravity square problem with quadratic elements
    T_pU_qall, _ = PUTQ(slip2dispstress, els, bcidxU, bcidxall, mu, nu)
    _, H_pT_qall = PUTQ(slip2dispstress, els, bcidxT, bcidxall, mu, nu)
    TH = [T_pU_qall ; alpha .* H_pT_qall] # Assemble combined linear operator

    # Particular solution and effective boundary conditions
    xnodes = transpose(els.xnodes[idx["B"], :])[:]
    ynodes = transpose(els.ynodes[idx["B"], :])[:]
    UB, _ = gravityparticularfunctions(xnodes, ynodes, g, rho, lambda, mu)    
    bcs = zeros(6 * els.endidx)
    bcs[1:2:6*nels] = UB[:, 1] # Bottom boundary (x-component)
    bcs[2:2:6*nels] = UB[:, 2] # Bottom boundary (y-component)
    bcs *= -1 # This is neccesary for the right answer and is consistent with derivation
    
    # BEM solve to get particular solution
    Ueffparticular = inv(TH) * bcs

    Uinteriorcomplementary, Sinteriorcomplementary = quaddispstress(slip2dispstress, x, y, els, bcidxall, quadstack(Ueffparticular[1:2:end]), quadstack(Ueffparticular[2:2:end]), mu, nu)
    Uinteriorparticular, Sinteriorparticular = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    U = @. Uinteriorcomplementary + Uinteriorparticular
    S = @. Sinteriorcomplementary + Sinteriorparticular
    Umag = sqrt.(U[:, 1].^2 + U[:, 2].^2)
    
    subplot(1, 2, 2)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts).+2*width, reshape(Umag, npts, npts), 50)
    colorbar()
    gca().set_xlim([-width, width])
    gca().set_ylim([0, 2*width])
    gca().set_aspect("equal")
    title("BEM solution")
end
femvsbemgravitybox()
