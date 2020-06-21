using Revise
using PyPlot
using Infiltrator
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
    # TODO: Quadratic elements
    # TODO: Move particular solution to Bem2d.jl
    # TODO: Move generation of modified BCs to function
    # TODO: Create a convenciece function for simple elements (addelsez?)
    # TODO: It's strange that the top of the model has to be at zero.  Can we generalize this?
    # TODO: Rule of thumb for choosing precondtioner value (alpha)?
    # TODO: This file should be renamed to gravitysquareparticular
    # TODO: Remove all other gravitysquareparticular functions

    close("all")
    alpha = 7e-8 # scalar preconditioner for traction terms
    fontsize = 20
    mu = 3e10
    lambda = 3e10
    nu = 0.25
    rho = 2700
    g = 9.81
    nels = 20
    npts = 25
    L = 1e4
    offset = 1000
    x, y = obsgrid(-L+offset, -2*L+offset, L-offset, 0-offset, npts) 

    # BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-L, -2*L, L, -2*L, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "B"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(L, -2*L, L, 0, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "R"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(L, 0, -L, 0, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "T"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(-L, 0, -L, -2*L, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "L"
    end
    standardize_elements!(els)

    # Common indexing
    idx = getidxdict(els)
    Bidx = idx["B"]
    RTLidx = [idx["R"] ; idx["T"]; idx["L"]]
    BRTLidx = collect(1:1:els.endidx)

    # Constant kernel matrices
    T_pB_qBRTL, H_pB_qBRTL = PUTC(slip2dispstress, els, Bidx, BRTLidx, mu, nu)
    T_pRTL_qBRTL, H_pRTL_qBRTL = PUTC(slip2dispstress, els, RTLidx, BRTLidx, mu, nu)

    # Trying a bottom BC only case
    bcseff = zeros(2 * els.endidx)
    UB, _ = gravityparticularfunctions(els.xcenter[idx["B"]], els.ycenter[idx["B"]], g, rho, lambda, mu)
    bcseff[@. 2*idx["B"]-1] = UB[:, 1] # Bottom boundary (x-component)
    bcseff[2*idx["B"]] = UB[:, 2] # Bottom boundary (y-component)
    bcseff[@. 2*idx["R"]-1] .= 0 # Right boundary (x-component)
    bcseff[2*idx["R"]] .= 0 # Right boundary (y-component)
    bcseff[@. 2*idx["T"]-1] .= 0 # Top boundary (x-component)
    bcseff[2*idx["T"]] .= 0 # Top boundary (y-component)    
    bcseff[@. 2*idx["L"]-1] .= 0 # Left boundary (x-component)
    bcseff[2*idx["L"]] .= 0 # Left boundary (y-component)
    bcseff *= -1 # This is neccesary for the right answer and is consistent with derivation

    # Solve BEM problem
    TH = [T_pB_qBRTL ; alpha .* H_pRTL_qBRTL]
    Ueffparticular = inv(TH) * bcseff
    
    # Quadratic kernel matrices
    T_pB_qBRTL_Q, H_pB_qBRTL_Q = PUTQ(slip2dispstress, els, Bidx, BRTLidx, mu, nu)
    T_pRTL_qBRTL_Q, H_pRTL_qBRTL_Q = PUTQ(slip2dispstress, els, RTLidx, BRTLidx, mu, nu)
    bcsQ = zeros(6 * els.endidx)
    idxQ = collect(1:1:length(bcsQ))    
    # bcsQ[2*idx["B"].-1] = @. -lambda*rho*g * els.xcenter[idx["B"]]*els.ycenter[idx["B"]] / (4*mu*(lambda+mu)) # Bottom boundary (x-component)
    # bcsQ[2*idx["B"]] = @. (rho*g) * (lambda*els.xcenter[idx["B"]]^2 + (lambda+2*mu)*els.ycenter[idx["B"]]^2) / (8*mu*(lambda+mu)) # Bottom boundary (y-component)
    bcsQ[idxQ[els.endidx + 1 : 2 : 2 * (els.endidx)]] .= 0 # Right boundary (x-component)
    bcsQ[@. idxQ[els.endidx + 1 : 2 : 2 * (els.endidx)] + 1] .= 0 # Right boundary (y-component)
    bcsQ[idxQ[2*els.endidx + 1 : 2 : 3 * (els.endidx)]] .= 0 # Top boundary (x-component)
    bcsQ[@. idxQ[2*els.endidx + 1 : 2 : 3 * (els.endidx)] + 1] .= 0 # Top boundary (y-component)    
    # bcsQ[2*idx["L"].-1] .= 0 # Left boundary (x-component)
    # bcsQ[2*idx["L"]] .= 0 # Left boundary (y-component)
    bcsQ *= -1 # This is neccesary for the right answer and is consistent with derivation
    TH = [T_pB_qBRTL_Q ; alpha .* H_pRTL_qBRTL_Q]
    UeffparticularQ = inv(TH) * bcsQ


    # Complementary solution from Ueff
    UinteriorBRTL, SinteriorBRTL = constdispstress(slip2dispstress, x, y, els, BRTLidx, Ueffparticular[1:2:end], Ueffparticular[2:2:end], mu, nu)
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), UinteriorBRTL, SinteriorBRTL, "ueff solution")

    # Particular solution
    Uinteriorparticular, Sinteriorparticular = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), Uinteriorparticular, Sinteriorparticular, "Particular solution")
    
    # Combine solutions following equations 14 in Pape and Bannerjee
    U = @. UinteriorBRTL + Uinteriorparticular
    S = @. SinteriorBRTL + Sinteriorparticular

    # Plot total solution particular + complementary
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), U, S, "uc + up")
end
gravitysquareparticular()
