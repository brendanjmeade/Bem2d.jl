using Revise
using PyPlot
using Bem2d


"""
    dislocationinabox()

Comparing half-space and dislocaiton in a box solutions
"""
function dislocationinabox()
    mu = 30e9
    nu = 0.25

    els = Elements(Int(1e5))
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-50, 0, 50, 0, nfreesurf) # Free surface
    addelsez!(els, x1, y1, x2, y2, "freesurf")
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-1, -1, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(els, x1, y1, x2, y2, "fault")

    # Conviencience structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)

    # Constant slip fault
    partialsconst["disp"]["fault"]["freesurf"], partialsconst["stress"]["fault"]["freesurf"], partialsconst["trac"]["fault"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    partialsconst["disp"]["freesurf"]["freesurf"], partialsconst["stress"]["freesurf"]["freesurf"], partialsconst["trac"]["freesurf"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    faultslipconst = sqrt(2) / 2 * [1 ; 1]
    ufullspaceconst = partialsconst["disp"]["fault"]["freesurf"] * faultslipconst
    ufreesurfaceconst = inv(partialsconst["trac"]["freesurf"]["freesurf"]) * (partialsconst["trac"]["fault"]["freesurf"] * faultslipconst)
    xplotconst = els.xcenter[idx["freesurf"]]


    fontsize = 20
    markersize = 12
    linewidth = 2.0
    close("all")
    figure(figsize = (12, 12))

    ax = subplot(2, 1, 1)
    plot(xplotconst, ufreesurfaceconst[1:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$u_x$ (m)", fontsize=fontsize)

    ax = subplot(2, 1, 2)
    plot(xplotconst, ufreesurfaceconst[2:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$u_y$ (m)", fontsize=fontsize)
    show()
end
dislocationinabox()
