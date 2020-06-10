using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d


"""
    gravitysquareparticular()

Experiments with gravity body force.
"""
function gravitysquareparticular()
    close("all")
    fontsize = 20
    mu = 3e10
    lambda = 3e10
    nu = 0.25
    rho = 2700
    g = 9.81
    nels = 20
    npts = 25
    L = 1e4
    x, y = obsgrid(-L+1, 0+1, L-1, 2*L-1, npts)

    #! BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "B"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(L, 0, L, 2*L, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "R"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(L, 2*L, -L, 2*L, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "T"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(-L, 2*L, -L, 0, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "L"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Kernels matrices (boundary -> boundary)
    Bidx = idx["B"]
    RTLidx = [idx["R"] ; idx["T"]; idx["L"]]
    BRTLidx = collect(1:1:els.endidx)
    T_pB_qBRTL, H_pB_qBRTL = PUTC(slip2dispstress, els, Bidx, BRTLidx, mu, nu)
    T_pRTL_qBRTL, H_pRTL_qBRTL = PUTC(slip2dispstress, els, RTLidx, BRTLidx, mu, nu)

    #! Particular solution technique with modified boundary conditoions
    alpha = 7e-8
    bcs = zeros(2 * 4 * nels)
    bcs[1:2:40] = @. -lambda*rho*g * els.xcenter[idx["B"]] * els.ycenter[idx["B"]] / (4*mu*(lambda + mu)) # Bottom boundary (x-component)
    bcs[2:2:40] = @. (rho*g) * (lambda*els.xcenter[idx["B"]]^2 + (lambda+2*mu)*els.ycenter[idx["B"]]^2) / (8*mu*(lambda + mu)) # Bottom boundary (y-component)
    bcs[41:2:80] .= 0 # Right boundary (x-component)
    bcs[42:2:80] = @. alpha * rho * g * (els.ycenter[idx["R"]]) # Right boundary (y-component)
    bcs[81:2:120] .= 0 # Top boundary (x-component)
    bcs[82:2:120] .= @. alpha * rho * g * (els.ycenter[idx["T"]]) # Top boundary (y-component)    
    bcs[121:2:160] .= 0 # Left boundary (x-component)
    bcs[122:2:160] = @. alpha * rho * g * (els.ycenter[idx["L"]])  # Left boundary (y-component)
    bcs *= -1 # Flip sign of all bcs?

    figure(figsize=(8, 8))
    subplot(2, 1, 1)
    plot(bcs[1:2:end])
    xlabel("idx")
    ylabel("x bcs")
    subplot(2, 1, 2)
    plot(bcs[2:2:end])
    xlabel("idx")
    ylabel("y bcs")

    TH = [T_pB_qBRTL ; alpha .* H_pRTL_qBRTL]
    Ueffparticular = inv(TH) * bcs
    
    #! Interior displacements from boundaries (Ueff)
    UinteriorBRTL, SinteriorBRTL = constdispstress(slip2dispstress, x, y, els, BRTLidx, Ueffparticular[1:2:end], Ueffparticular[2:2:end], mu, nu)
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), UinteriorBRTL, SinteriorBRTL, "Particular integral method")
    
    #! Single quiver plot of interior displacements
    figure(figsize=(8, 8))
    quiver(x, y, UinteriorBRTL[:, 1], UinteriorBRTL[:, 2])
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    title("internal displacements", fontsize=fontsize)
end
gravitysquareparticular()
