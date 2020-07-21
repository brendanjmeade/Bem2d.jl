using Revise
using PyPlot
using Infiltrator
using Bem2d


"""
    dislocationinabox()

Comparing half-space and dislocaiton in a box solutions
"""
function dislocationinabox()
    close("all")
    mu = 30e9
    nu = 0.25

    # Elment geometries and data structures for half space approximation
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-50e3, 0, 50e3, 0, 200) # Free surface
    addelsez!(els, x1, y1, x2, y2, "surf")
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, 1) # 45 degree dipping fault
    addelsez!(els, x1, y1, x2, y2, "fault")
    idx = getidxdict(els)

    # Solve the BEM problem
    _, H_surf_fault = PUTC(slip2dispstress, els, idx["surf"], idx["fault"], mu, nu)
    _, H_surf_surf = PUTC(slip2dispstress, els, idx["surf"], idx["surf"], mu, nu)
    faultslip = sqrt(2) / 2 * [1 ; 1]
    Ueff = inv(H_surf_surf) * (H_surf_fault * faultslip)
    
    # Evaluate in a volume
    npts = 50
    x, y = obsgrid(-50e3, -20e3, 50e3, 0, npts)
    Ufault, Sfault = constdispstress(slip2dispstress, x, y, els, idx["fault"], 1, 1, mu, nu)
    Usurf, Ssurf = constdispstress(slip2dispstress, x, y, els, idx["surf"], Ueff[1:2:end], Ueff[2:2:end], mu, nu)
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), Ufault, Sfault, "fault only")
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), Usurf, Ssurf, "surface response")
    
    # Plot surface displacements
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    figure(figsize = (12, 12))
    ax = subplot(2, 1, 1)
    plot(els.xcenter[idx["surf"]], Ueff[1:2:end], "b.", label="halfspace")
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([-1.0, 1.0])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$u_x$ (m)", fontsize=fontsize)

    ax = subplot(2, 1, 2)
    plot(els.xcenter[idx["surf"]], Ueff[2:2:end], "b.", label="halfspace")
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([-1.0, 1.0])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$u_y$ (m)", fontsize=fontsize)
    show()
end
dislocationinabox()
