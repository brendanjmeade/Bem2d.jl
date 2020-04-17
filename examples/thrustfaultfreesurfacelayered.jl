using Revise
using PyCall
using PyPlot
using Bem2d
ow = pyimport("okada_wrapper") # from okada_wrapper import dc3dwrapper


"""
    okadaexample()

Reference Okada calculation for a thrust fault tipping at at 45 degrees.
"""
function okadaexample()
    xokada = collect(LinRange(-5, 5, 1000))
    uxokada = zeros(length(xokada))
    uyokada = zeros(length(xokada))

    for i in 1:length(xokada)
        # Fault dipping at 45 degrees
        _, u, _ = ow.dc3dwrapper(
            2.0 / 3.0,
            [0, xokada[i] + 0.5, 0],
            0.5,
            45,  # 135
            [-1000, 1000],
            [-sqrt(2) / 2, sqrt(2) / 2],
            [0.0, 1.0, 0.0],
        )
        uxokada[i] = u[2]
        uyokada[i] = u[3]
    end
    return xokada, uxokada, uyokada
end


"""
    thrustfaultfreesurfacelayered()

Comparison of surface displacements near a thrust fault dipping
at 45 degrees.  Includes both constant and quadratic elements.
"""
function thrustfaultfreesurfacelayered()
    mu = 30e9
    nu = 0.25

    #! Free surface
    els = Elements(Int(1e5))
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-10, 0, 10, 0, nfreesurf)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "surface"
    end
    standardize_elements!(els)

    #! Buried interface
    ninterface = 200
    x1, y1, x2, y2 = discretizedline(-10, -2, 10, -2, ninterface)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "interface"
    end
    standardize_elements!(els)

    #! 45 degree dipping fault
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-1, -1, 0, 0, nfault)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Kernels
    Tfaultsurface, _,
    Hfaultsurface = partialsconstdispstress(slip2dispstress, els,
                                            idx["fault"], idx["surface"],
                                            mu, nu)

    Tsurfacesurface, _,
    Hsurfacesurface = partialsconstdispstress(slip2dispstress, els,
                                              idx["surface"], idx["surface"],
                                              mu, nu)

    Tfaultinterface, _,
    Hfaultinterface = partialsconstdispstress(slip2dispstress, els,
                                              idx["fault"], idx["interface"],
                                              mu, nu)
    
    Tinterfaceinterface, _,
    Hinterfaceinterface = partialsconstdispstress(slip2dispstress, els,
                                                  idx["interface"], idx["interface"],
                                                  mu, nu)
    

    #! Solve the BEM problem
    faultslip = sqrt(2) / 2 * [1 ; 1]
    ufreesurface = inv(Hsurfacesurface) * (Hfaultsurface * faultslip)
    xplot = els.xcenter[idx["surface"]]

    #! Stress and displacement continuity at buried interface

    #! Okada halfspace example for reference
    xokada, uxokada, uyokada = okadaexample()

    #! Plot surface displacement profiles
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    close("all")
    figure(figsize = (12, 12))

    ax = subplot(2, 1, 1)
    plot(xokada, uxokada, "-k", linewidth=linewidth, label="Okada")
    plot(xplot, ufreesurface[1:2:end], "r.", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([-1.00, 0.0, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize=fontsize)
    ylabel(L"$u_x$ (m)", fontsize=fontsize)

    ax = subplot(2, 1, 2)
    plot(xokada, uyokada, "-k", linewidth=linewidth, label="Okada")
    plot(xplot, ufreesurface[2:2:end], "r.", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([-1.00, 0.0, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize=fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$u_y$ (m)", fontsize=fontsize)
    show()
end
thrustfaultfreesurfacelayered()
