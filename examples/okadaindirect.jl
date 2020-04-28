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
    thrustfaultfreesurface()

Comparison of surface displacements near a thrust fault dipping
at 45 degrees.  Includes both constant and quadratic elements.
"""
function thrustfaultfreesurface()
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
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-50, 0, 50, 0, nfreesurf)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "surface"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Hacky old version
    Tfaultfreesurf, _, Hfaultfreesurf = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["surface"], mu, nu)
    Tfreesurffreesurf, _, Hfreesurffreesurf = partialsconstdispstress(slip2dispstress, els, idx["surface"], idx["surface"], mu, nu)
    faultslip = sqrt(2) / 2 * [1 ; 1]
    ufullspaceconst = Tfaultfreesurf * faultslip
    ufreesurfaceconst = inv(Hfreesurffreesurf) * (Hfaultfreesurf * faultslip)
    xplotconst = els.xcenter[idx["surface"]]

    # #! Formal indirect version
    Tfaultfault, _, _ = PUSTC(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    _, _, Hfaultfreesurf = PUSTC(slip2dispstress, els, idx["fault"], idx["surface"], mu, nu)
    Tfreesurffault, _, _ = PUSTC(slip2dispstress, els, idx["surface"], idx["fault"], mu, nu)
    Tfreesurffreesurf, _, Hfreesurffreesurf = PUSTC(slip2dispstress, els, idx["surface"], idx["surface"], mu, nu)
    mat = zeros(2*els.endidx, 2*els.endidx)
    mat[1:2, 1:2] = Tfaultfault
    mat[1:2, 3:end] = Hfaultfreesurf
    mat[3:end, 1:2] = Tfreesurffault
    mat[3:end, 3:end] = Hfreesurffreesurf
    bcs = zeros(2*els.endidx)
    bcs[1:2] = faultslip
    ueff = inv(mat) * bcs
    matold = mat

    #! Formal indirect version
    # Tfault, _, Hfault = PUSTC(slip2dispstress, els, idx["fault"], collect(1:1:els.endidx), mu, nu)
    # Tsurface, _, Hsurface = PUSTC(slip2dispstress, els, idx["freesurf"], collect(1:1:els.endidx), mu, nu)

    # @show size(Tfault)
    # @show size(Hfault)
    # @show size(Tsurface)
    # @show size(Hsurface)

    # mat = zeros(2*els.endidx, 2*els.endidx)
    # mat[1:2, :] = Tfault
    # mat[3:end, :] = Hsurface
    # bcs = zeros(2*els.endidx)
    # bcs[1:2] = faultslipconst
    # ueff = inv(mat) * bcs

    # @infiltrate
    # return

    # Forward evaluation at free surface
    matinterior = zeros(2*length(idx["surface"]), 2*els.endidx)
    matinterior[:, 1:2] = Tfreesurffault
    matinterior[:, 3:end] = Tfreesurffreesurf
    usurf = matinterior * ueff

    #! Okada solution
    xokada = collect(LinRange(-5, 5, 1000))
    uxokada, uyokada = okadalocal(xokada)

    #! Plot comparison between BEM and analytic
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    close("all")
    figure(figsize = (12, 12))

    ax = subplot(2, 1, 1)
    plot(xokada, uxokada, "-k", linewidth=linewidth, label="Okada")
    plot(xplotconst, ufreesurfaceconst[1:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label = "hacky BEM")
    plot(xplotconst, usurf[1:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "indirect BEM")
    ylabel(L"$u_x$ (m)", fontsize=fontsize)
    plotformat(fontsize)

    ax = subplot(2, 1, 2)
    plot(xokada, uyokada, "-k", linewidth=linewidth, label="Okada")
    plot(xplotconst, ufreesurfaceconst[2:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label = "hacky BEM")
    plot(xplotconst, usurf[2:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "indirect BEM")
    ylabel(L"$u_y$ (m)", fontsize=fontsize)
    plotformat(fontsize)
    show()
end
thrustfaultfreesurface()
