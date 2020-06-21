using Revise
using PyPlot
using Infiltrator
using LinearAlgebra
using Bem2d


"""
    gravityparticularfunctions()

From Pape and Bannerjee 1987
"""
function gravityparticularfunctions(x, y, g, rho, lambda, mu)
    U = zeros(length(x), 2)
    S = zeros(length(x), 3)
    U[:, 1] = @. -lambda * rho * g / (4 * mu * (lambda + mu)) * x * y  # Pape and Banerjee (1987) equation (6a)
    U[:, 2] = @. (rho * g) / (8 * mu * (lambda + mu)) * (lambda * x^2 + (lambda + 2 * mu) * y^2) # Pape and Banerjee (1987) equation (6b)
    S[:, 1] .= 0 # Pape and Banerjee (1987) equation (6c)
    S[:, 2] = @. rho * g * y  # Pape and Banerjee (1987) equation (6d)
    S[:, 3] .= 0 # Pape and Banerjee (1987) equation (6e)
    return U, S
end


"""
    gravitysquareparticular()

Experiments with gravity body force.
"""
function gravitysquareparticular()
    # TODO: Move particular solution to Bem2d.jl
    # TODO: It's strange that the top of the model has to be at zero.  Can we generalize this?
    # TODO: Rule of thumb for choosing precondtioner value (alpha)?
    # TODO: Make a version of gravityparticularfunctions() for BC generation

    close("all")
    alpha = 7e-8 # scalar preconditioner for traction terms
    fontsize = 20
    mu = 3e10
    lambda = 3e10
    nu = 0.25
    rho = 2700
    g = 9.81
    nels = 20
    npts = 100
    L = 1e4
    offset = 10
    x, y = obsgrid(-L+offset, -2*L+offset, L-offset, 0-offset, npts) 

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-L, -2*L, L, -2*L, nels) # Bottom
    addelsez!(els, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(L, -2*L, L, 0, nels) # Right hand side
    addelsez!(els, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(L, 0, -L, 0, nels) # Top
    addelsez!(els, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-L, 0, -L, -2*L, nels) # Left hand side
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

    # Evaluate and plot interior solution
    Uinteriorcomplementary, Sinteriorcomplementary = quaddispstress(slip2dispstress, x, y, els, bcidxall, quadstack(Ueffparticular[1:2:end]), quadstack(Ueffparticular[2:2:end]), mu, nu)
    Uinteriorparticular, Sinteriorparticular = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    U = @. Uinteriorcomplementary + Uinteriorparticular
    S = @. Sinteriorcomplementary + Sinteriorparticular
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), Uinteriorcomplementary, Sinteriorcomplementary, "Complementary solution")
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), Uinteriorparticular, Sinteriorparticular, "Particular solution")
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), U, S, "Complementary + Particular solutions")
end
gravitysquareparticular()
