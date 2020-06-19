using Revise
using PyPlot
using Bem2d


"""
    gravitysquareparticular()

Experiments with gravity body force.
"""
function gravitysquareparticular()
    close("all")
    alpha = 7e-8 # scalar precondtioner for traction terms
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
    x, y = obsgrid(-L+offset, 0+offset, L-offset, 2*L-offset, npts)

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

    # Particular solution technique with modified boundary conditoions    
    # bcseff = zeros(2 * 4 * nels)
    # bcseff[1:2:40] = @. -lambda*rho*g * els.xcenter[idx["B"]]*els.ycenter[idx["B"]] / (4*mu*(lambda + mu)) # Bottom boundary (x-component)
    # bcseff[2:2:40] = @. (rho*g) * (lambda*els.xcenter[idx["B"]]^2 + (lambda+2*mu)*els.ycenter[idx["B"]]^2) / (8*mu*(lambda + mu)) # Bottom boundary (y-component)
    # bcseff[41:2:80] .= 0 # Right boundary (x-component)
    # bcseff[42:2:80] .= 0 # Right boundary (y-component)
    # bcseff[81:2:120] .= 0 # Top boundary (x-component)
    # bcseff[82:2:120] = @. alpha * rho * g * (els.ycenter[idx["T"]]) # Top boundary (y-component)    
    # bcseff[121:2:160] .= 0 # Left boundary (x-component)
    # bcseff[122:2:160] .= 0 # Left boundary (y-component)
    # bcseff *= -1 # This is neccesary for the right answer and is consistent with derivation

    # Tryin a bottom BC only case
    bcseff = zeros(2 * 4 * nels)
    # bcseff[1:2:40] = @. -lambda*rho*g * els.xcenter[idx["B"]]*els.ycenter[idx["B"]] / (4*mu*(lambda + mu)) # Bottom boundary (x-component)
    # bcseff[2:2:40] = @. (rho*g) * (lambda*els.xcenter[idx["B"]]^2 + (lambda+2*mu)*els.ycenter[idx["B"]]^2) / (8*mu*(lambda + mu)) # Bottom boundary (y-component)
    bcseff[1:2:40] = @. -lambda*rho*g * els.xcenter[idx["B"]]*-20000 / (4*mu*(lambda+mu)) # Bottom boundary (x-component)
    bcseff[2:2:40] = @. (rho*g) * (lambda*els.xcenter[idx["B"]]^2 + (lambda+2*mu)*20000^2) / (8*mu*(lambda+mu)) # Bottom boundary (y-component)
    bcseff[41:2:80] .= 0 # Right boundary (x-component)
    bcseff[42:2:80] .= 0 # Right boundary (y-component)
    bcseff[81:2:120] .= 0 # Top boundary (x-component)
    # bcseff[82:2:120] = @. alpha * rho * g * (els.ycenter[idx["T"]]) # Top boundary (y-component)    
    bcseff[82:2:120] .= 0 # Top boundary (y-component)    
    bcseff[121:2:160] .= 0 # Left boundary (x-component)
    bcseff[122:2:160] .= 0 # Left boundary (y-component)
    bcseff *= -1 # This is neccesary for the right answer and is consistent with derivation



    # Solve BEM problem
    TH = [T_pB_qBRTL ; alpha .* H_pRTL_qBRTL]
    Ueffparticular = inv(TH) * bcseff
    
    # Interior displacements from Ueff
    UinteriorBRTL, SinteriorBRTL = constdispstress(slip2dispstress, x, y, els, BRTLidx, Ueffparticular[1:2:end], Ueffparticular[2:2:end], mu, nu)
    
    # Plot 6 panel interior evaluation
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), UinteriorBRTL, SinteriorBRTL, "Particular integral method")

    # Single quiver plot of interior displacements
    figure(figsize=(16, 8))
    subplot(2, 2, 1)
    plot(bcseff[1:2:end], "-r", label="bcs")
    plot(Ueffparticular[1:2:end], "-k", label="ueff")
    text(10, 0.0, "bottom", horizontalalignment="center", verticalalignment="center")
    text(30, 0.0, "right", horizontalalignment="center", verticalalignment="center")
    text(50, 0.0, "top", horizontalalignment="center", verticalalignment="center")
    text(70, 0.0, "left", horizontalalignment="center", verticalalignment="center")
    legend()
    xlabel(" ")
    ylabel("value")
    title("x boundary conditions")

    subplot(2, 2, 3)
    plot(bcseff[2:2:end], "-r", label="bcs")
    plot(Ueffparticular[2:2:end], "-k", label="ueff")
    text(10, 0.0, "bottom", horizontalalignment="center", verticalalignment="center")
    text(30, 0.0, "right", horizontalalignment="center", verticalalignment="center")
    text(50, 0.0, "top", horizontalalignment="center", verticalalignment="center")
    text(70, 0.0, "left", horizontalalignment="center", verticalalignment="center")
    legend()
    xlabel("index")
    ylabel("value")
    title("y boundary conditions")

    subplot(1, 2, 2)
    quiver(x, y, UinteriorBRTL[:, 1], UinteriorBRTL[:, 2])
    xlabel("x (m)")
    ylabel("y (m)")
    gca().set_aspect("equal")
    title("internal displacements from ueff")
end
gravitysquareparticular()
