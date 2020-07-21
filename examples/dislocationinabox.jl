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
    x1, y1, x2, y2 = discretizedline(-50e3, 0, 50e3, 0, 50) # Free surface
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
    Ufault, _ = constdispstress(slip2dispstress, els.xcenter[idx["surf"]], els.ycenter[idx["surf"]], els, idx["fault"], 1, 1, mu, nu)
    Usurf, _ = constdispstress(slip2dispstress, els.xcenter[idx["surf"]], els.ycenter[idx["surf"]], els, idx["surf"], Ueff[1:2:end], Ueff[2:2:end], mu, nu)
    Utotal = @. Ufault - Usurf
    
    # Plot surface displacements
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    figure(figsize = (12, 12))
    ax = subplot(2, 1, 1)
    plot(els.xcenter[idx["surf"]], Ueff[1:2:end], "b.", label="Ueff")
    plot(els.xcenter[idx["surf"]], Ufault[:, 1], "r+", label="fault only")
    plot(els.xcenter[idx["surf"]], Usurf[:, 1], "rx", label="fault only")
    plot(els.xcenter[idx["surf"]], Utotal[:, 1], "r.", label="total")
    
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
