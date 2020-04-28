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
    mu = 30e9
    nu = 0.25

    # 45 degree dipping fault
    els = Elements(Int(1e5))
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
    faultslip = sqrt(2) / 2 * [1 ; 1]
    xbem = els.xcenter[idx["surface"]]

    #! Hacky BEM
    ufullspacehacky, uhalfspacehacky = hackybem(els, idx, faultslip, mu, nu)

    #! Formal indirect BEM
    T_pfault_qfault, _, _ = PUSTC(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    _, _, H_pfault_qsurface = PUSTC(slip2dispstress, els, idx["fault"], idx["surface"], mu, nu)
    T_psurface_qfault, _, _ = PUSTC(slip2dispstress, els, idx["surface"], idx["fault"], mu, nu)
    T_psurface_qsurface, _, H_psurface_qsurface = PUSTC(slip2dispstress, els, idx["surface"], idx["surface"], mu, nu)
    mat = zeros(2*els.endidx, 2*els.endidx)
    mat[1:2, 1:2] = T_pfault_qfault
    mat[3:end, 1:2] = T_psurface_qfault #TODO: I think this should be T_pfault_qsurface
    mat[1:2, 3:end] = H_pfault_qsurface  #TODO: I think this should be T_psurface_qfault
    mat[3:end, 3:end] = H_psurface_qsurface
    bcs = zeros(2*els.endidx)
    bcs[1:2] = faultslip
    ueff = inv(mat) * bcs

    # Interior evaluation at free surface
    matinterior = zeros(2*length(idx["surface"]), 2*els.endidx)
    matinterior[:, 1:2] = T_psurface_qfault
    matinterior[:, 3:end] = T_psurface_qsurface
    @show size(T_psurface_qfault)
    ufullspaceindirect = T_psurface_qfault * (ueff[1:2])
    uhalfspaceindirect = matinterior * ueff

    #! Okada solution
    xokada = collect(LinRange(-5, 5, 1000))
    uxokada, uyokada = okadalocal(xokada)

    #! Plot comparison between BEM and Okada
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    close("all")
    figure(figsize = (20, 20))

    ax = subplot(2, 1, 1)
    plot(xokada, uxokada, "-k", linewidth=linewidth, label="Okada (halfspace)")
    plot(xbem, ufullspacehacky[1:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label="hacky BEM (fullspace)")
    plot(xbem, uhalfspacehacky[1:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label="hacky BEM (halfspace)")
    plot(xbem, ufullspaceindirect[1:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label="indirect BEM (fullspace)")
    plot(xbem, uhalfspaceindirect[1:2:end], "rx", markeredgewidth=linewidth, markersize=markersize, label="indirect BEM (halfspace)")
    ylabel(L"$u_x$ (m)", fontsize=fontsize)
    plotformat(fontsize)

    ax = subplot(2, 1, 2)
    plot(xokada, uyokada, "-k", linewidth=linewidth, label="Okada (halfspace)")
    plot(xbem, ufullspacehacky[2:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label = "hacky BEM (fullspace)")
    plot(xbem, uhalfspacehacky[2:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label = "hacky BEM (halfspace)")
    plot(xbem, ufullspaceindirect[2:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label="indirect BEM (fullspace)")
    plot(xbem, uhalfspaceindirect[2:2:end], "rx", markeredgewidth=linewidth, markersize=markersize, label = "indirect BEM (halfspace)")
    ylabel(L"$u_y$ (m)", fontsize=fontsize)
    plotformat(fontsize)
    show()

    @infiltrate
end
okadaindirect()
