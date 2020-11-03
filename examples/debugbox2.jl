using Revise
using PyPlot
using Infiltrator
using LinearAlgebra
using Bem2d


"""
    bemsolve(g, rho, lambda, mu)

Solve BEM with particular integral approach with gravitysquareparticular
"""
function bemsolve(els, nels, g, rho, lambda, mu, nu, x, y)
    # Build BEM operator, TH
    idx = getidxdict(els)
    T_pU_qall, _ = PUTC(slip2dispstress, els, idx["B"], 1:1:els.endidx, mu, nu)
    _, H_pT_qall = PUTC(slip2dispstress, els, [idx["R"] ; idx["T"]; idx["L"]], 1:1:els.endidx, mu, nu)
    TH = [T_pU_qall ; H_pT_qall]

    # Modified displacment BCs for bottom
    Upart, _ = gravityparticularfunctions(els.xcenter[idx["B"]], els.ycenter[idx["B"]], g, rho, lambda, mu)
    bcs = zeros(2 * els.endidx)
    bcs[1:2:2*nels] = Upart[:, 1]
    bcs[2:2:2*nels] = Upart[:, 2]

    # Modified traction BCs for top
    bcs[4*nels+2:2:6*nels] = rho * g * els.ynormal[idx["T"]]
    bcs *= -1
    Ueff = TH \ bcs
    
    # Forward solution on grid
    Ucomp, _ = constdispstress(slip2dispstress, x, y, els, [idx["B"] ; idx["R"] ; idx["T"] ; idx["L"]],
                               Ueff[1:2:end], Ueff[2:2:end], mu, nu)
    Uint, _ = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    Utotal = Ucomp .+ Uint
    Umag = sqrt.(Utotal[:, 1].^2 + Utotal[:, 2].^2)
    return Umag
end


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
    npts = 100
    L = 1e4
    offset = 10

    ##### Below ground box
    x, y = obsgrid(-L+offset, -2*L+offset, L-offset, 0-offset, npts) 
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-L, -2*L, L, -2*L, nels) # Bottom
    addelsez!(els, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(L, -2*L, L, 0, nels) # Right hand side
    addelsez!(els, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(L, 0, -L, 0, nels) # Top
    addelsez!(els, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-L, 0, -L, -2*L, nels) # Left hand side
    addelsez!(els, x1, y1, x2, y2, "L")

    Umag = bemsolve(els, nels, g, rho, lambda, mu, nu, x, y)
    figure(figsize=(15, 5))
    subplot(1, 3, 1)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Umag, npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("y(top) = 0")

    ##### Now with top of box not at zero
    xoffset = 0 # Horizontal box translation...works fine
    yoffset = 2 * L
    x .+= xoffset
    els.x1 .+= xoffset
    els.x2 .+= xoffset
    els.xcenter .+= xoffset
    y .+= yoffset
    els.y1 .+= yoffset
    els.y2 .+= yoffset
    els.ycenter .+= yoffset

    Umagoffset = bemsolve(els, nels, g, rho, lambda, mu, nu, x, y)
    subplot(1, 3, 2)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Umagoffset, npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("y(top) = yoffset")

    subplot(1, 3, 3)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(abs.(Umagoffset - Umag), npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("residual")

end
gravitysquareparticular()
