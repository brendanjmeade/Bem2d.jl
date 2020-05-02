using Revise
using PyCall
using PyPlot
using Infiltrator
using Bem2d
ow = pyimport("okada_wrapper") # from okada_wrapper import dc3dwrapper


"""
    plotformat(fontsize)

Standard plot formatting.
"""
function plotformat(fontsize)
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([-1.00, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00])
    legend(fontsize=fontsize)
    gca().tick_params("both", labelsize=fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize)
    return nothing
end


"""
    okadalocal(xokada)

Okada profile for a 45 dipping fault.
"""
function okadalocal(xokada)
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
    return uxokada, uyokada
end


"""
    hackybem(els, mu, nu)

Hacky BEM approximation to Okada profile for a 45 dipping fault.
"""
function hackybem(els, idx, faultslip, mu, nu)
    Tfaultsurface, _, Hfaultsurface = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["surface"], mu, nu)
    Tsurfacesurface, _, Hsurfacesurface = partialsconstdispstress(slip2dispstress, els, idx["surface"], idx["surface"], mu, nu)
    ufullspacehacky = Tfaultsurface * faultslip
    usurfacehacky = inv(Hsurfacesurface) * (Hfaultsurface * faultslip)
    return ufullspacehacky, usurfacehacky
end


"""
    okadaindirect()

Comparison of surface displacements near a thrust fault dipping
at 45 degrees.  Includes both constant and quadratic elements.
"""
function okadaindirect()
    close("all")
    mu = 30e9
    nu = 0.25

    # 45 degree dipping fault
    els = Elements(Int(1e5))
    nfault = 10
    x1, y1, x2, y2 = discretizedline(-1, -1, 0, 0, nfault)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Free surface
    nfreesurf = 20
    x1, y1, x2, y2 = discretizedline(-5, 0, 5, 0, nfreesurf)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "surface"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Parameters for BEM solutions
    xbem = els.xcenter[idx["surface"]]

    #! Hacky BEM
    ufullspacehacky, uhalfspacehacky = hackybem(els, idx, sqrt(2)/2 .* ones(2*nfault), mu, nu)

    #! Formal indirect BEM
    _, _, H_psurface_qfault = PUSTC(slip2dispstress, els, idx["surface"], idx["fault"], mu, nu)
    _, _, H_psurface_qsurface = PUSTC(slip2dispstress, els, idx["surface"], idx["surface"], mu, nu)
    bcs = sqrt(2.0)/2.0 * ones(2*nfault)
    ueff = inv(H_psurface_qsurface) * H_psurface_qfault * bcs
    uhalfspaceindirect = ueff

    #! Okada solution
    xokada = collect(LinRange(-5, 5, 1000))
    uxokada, uyokada = okadalocal(xokada)

    #! Plot comparison between BEM and Okada
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    figure(figsize = (20, 20))

    ax = subplot(2, 1, 1)
    plot(xokada, uxokada, "-k", linewidth=linewidth, label="Okada (halfspace)")
    # plot(xbem, ufullspacehacky[1:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label="hacky BEM (fullspace)")
    plot(xbem, uhalfspacehacky[1:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label="hacky BEM (halfspace)")
    plot(xbem, uhalfspaceindirect[1:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label="indirect BEM (halfspace)")
    ylabel(L"$u_x$ (m)", fontsize=fontsize)
    plotformat(fontsize)

    ax = subplot(2, 1, 2)
    plot(xokada, uyokada, "-k", linewidth=linewidth, label="Okada (halfspace)")
    # plot(xbem, ufullspacehacky[2:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label = "hacky BEM (fullspace)")
    plot(xbem, uhalfspacehacky[2:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label = "hacky BEM (halfspace)")
    plot(xbem, uhalfspaceindirect[2:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "indirect BEM (halfspace)")
    ylabel(L"$u_y$ (m)", fontsize=fontsize)
    plotformat(fontsize)
    show()
    @infiltrate
end
okadaindirect()
