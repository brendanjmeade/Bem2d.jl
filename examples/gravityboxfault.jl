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
    gravitydislocation()

Experiments with gravity body force and a dislocation
"""
function gravitydislocation()
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
    offset = 0.1
    x, y = obsgrid(-30000+offset, -20000+offset, 30000-offset, 0-offset, npts) 

    # Define external boundary geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-30000, -20000, 30000, -20000, nels) # Bottom
    addelsez!(els, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(30000, -20000, 30000, 0, nels) # Right hand side
    addelsez!(els, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(30000, 0, -30000, 0, nels) # Top
    addelsez!(els, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-30000, 0, -30000, -20000, nels) # Left hand side
    addelsez!(els, x1, y1, x2, y2, "L")
    
    # Define dislocation geometry
    x1, y1, x2, y2 = discretizedline(-5000, -5000, 0, 0, 1)
    addelsez!(els, x1, y1, x2, y2, "D")

    # Common indexing
    idx = getidxdict(els)
    bcidxU = idx["B"] # Boundaries with *displacement* BCs
    bcidxT = [idx["R"] ; idx["T"]; idx["L"]] # Boundaries with *traction* BCs
    bcidxall = [idx["B"] ; idx["R"] ; idx["T"]; idx["L"]] # All exterior boundaries
    nbcels = length(bcidxall)

    ### Gravity only solution ###
    # Kernels
    T_pU_qall, H_pU_qall = PUTQ(slip2dispstress, els, bcidxU, bcidxall, mu, nu)
    T_pT_qall, H_pT_qall = PUTQ(slip2dispstress, els, bcidxT, bcidxall, mu, nu)
    THgravity = [T_pU_qall ; alpha .* H_pT_qall] # Assemble combined linear operator for gravity problem

    # Particular solution and effective boundary conditions
    xnodes = transpose(els.xnodes[idx["B"], :])[:]
    ynodes = transpose(els.ynodes[idx["B"], :])[:]
    UB, _ = gravityparticularfunctions(xnodes, ynodes, g, rho, lambda, mu)    
    bcsgravity = zeros(6 * nbcels)
    bcsgravity[1:2:6*nels] = UB[:, 1] # Bottom boundary (x-component)
    bcsgravity[2:2:6*nels] = UB[:, 2] # Bottom boundary (y-component)
    bcsgravity *= -1 # This is neccesary for the right answer and is consistent with derivation

    # BEM solve to get particular solution
    Ueffparticular = inv(THgravity) * bcsgravity

    # Evaluate and plot interior solution
    Uinteriorcomplementary, Sinteriorcomplementary = quaddispstress(slip2dispstress, x, y, els, bcidxall, quadstack(Ueffparticular[1:2:end]), quadstack(Ueffparticular[2:2:end]), mu, nu)
    Uinteriorparticular, Sinteriorparticular = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    U = @. Uinteriorcomplementary + Uinteriorparticular
    S = @. Sinteriorcomplementary + Sinteriorparticular
    # plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), U, S, "Complementary + Particular solutions")

    ### Dislocation only solution ###
    # Kernels
    T_pB_qD, H_pB_qD = PUTQ(slip2dispstress, els, idx["B"], idx["D"], mu, nu)
    T_pRTL_qD, H_pRTL_qD = PUTQ(slip2dispstress, els, bcidxT, idx["D"], mu, nu) 
    T_pB_qBRTL, H_pD_qBRTL = PUTQ(slip2dispstress, els, idx["B"], bcidxall, mu, nu)
    T_pRTL_qBRTL, H_pRTL_qBRTL = PUTQ(slip2dispstress, els, bcidxT, bcidxall, mu, nu)
    Ueffdislocation = -inv([T_pB_qBRTL ; alpha .* H_pRTL_qBRTL]) * [T_pB_qD ; alpha .* H_pRTL_qD]

    # Forward model for volume
    Uinteriorbcs, Sinteriorbcs = quaddispstress(slip2dispstress, x, y, els, bcidxall, quadstack(Ueffdislocation[1:2:end]), quadstack(Ueffdislocation[2:2:end]), mu, nu)
    Uinteriordislocation, Sinteriordislocation = quaddispstress(slip2dispstress, x, y, els, idx["D"], [1 1 1], [1 1 1], mu, nu)
    U = Uinteriordislocation .+ Uinteriorbcs
    S = Sinteriordislocation .+ Sinteriorbcs

    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts),
               Uinteriordislocation, Sinteriordislocation, "Dislocation only")
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts),
               Uinteriorbcs, Sinteriorbcs, "Boundary correction")
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts),
               U, S, "Dislocation + boundary correction")
end
gravitydislocation()
