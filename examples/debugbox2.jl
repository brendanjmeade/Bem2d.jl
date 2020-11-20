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
    bcs[1:2:2*nels] = -Upart[:, 1]
    bcs[2:2:2*nels] = -Upart[:, 2]

    # Modified traction BCs for top
    bcs[4*nels+2:2:6*nels] = -rho * g * (els.ycenter[idx["T"]] .* els.ynormal[idx["T"]])
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
    nels = 100
    npts = 200
    L = 1e4
    offset = 100

    ##### Below ground box
    top = 0
    x, y = obsgrid(-L+offset, top-2*L+offset, L-offset, top-offset, npts) 
    els = Elements(Int(1e5))
    addelsez!(els, discretizedline(-L, top-2*L, L, top-2*L, nels)..., "B") # Bottom
    addelsez!(els, discretizedline(L, top-2*L, L, top, nels)... , "R") # Right hand side
    addelsez!(els, discretizedline(L, top, -L, top, nels)..., "T") # Top
    addelsez!(els, discretizedline(-L, top, -L, top-2*L, nels)..., "L") # Left hand side

    Umag = bemsolve(els, nels, g, rho, lambda, mu, nu, x, y)
    figure(figsize=(15, 5))
    subplot(1, 3, 1)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Umag, npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("y(top) = 0")

    ##### Now with top of box not at zero
    top = 2*L
    x, y = obsgrid(-L+offset, top-2*L+offset, L-offset, top-offset, npts) 
    els = Elements(Int(1e5))
    addelsez!(els, discretizedline(-L, top-2*L, L, top-2*L, nels)..., "B") # Bottom
    addelsez!(els, discretizedline(L, top-2*L, L, top, nels)... , "R") # Right hand side
    addelsez!(els, discretizedline(L, top, -L, top, nels)..., "T") # Top
    addelsez!(els, discretizedline(-L, top, -L, top-2*L, nels)..., "L") # Left hand side

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
