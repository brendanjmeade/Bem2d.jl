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
    okadalocalinterior(xokada)

Okada profile for a 45 dipping fault.
"""
function okadalocalinterior(x, y)
    uxokada = zeros(length(x))
    uyokada = zeros(length(x))

    for i in 1:length(x)
        # Fault dipping at 45 degrees
        _, u, _ = ow.dc3dwrapper(
            2.0 / 3.0,
            [0, x[i] + 0.5, y[i]],
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
The purpose of keeping this is that while a little "hacky" it 
demonstrates the conceptual behind a mechanical interpretation of
the indirect BEM for the okada problem
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
    nfault = 20
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
    nfreesurf = 60
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
    T_psurface_qfault, H_psurface_qfault = PUTC(slip2dispstress, els, idx["surface"], idx["fault"], mu, nu)
    T_psurface_qsurface, H_psurface_qsurface = PUTC(slip2dispstress, els, idx["surface"], idx["surface"], mu, nu)
    bcs = sqrt(2.0)/2.0 * ones(2*nfault)
    ueff = inv(H_psurface_qsurface) * H_psurface_qfault * bcs
    uhalfspaceindirect = ueff

    #! Interior evaluation
    #TODO: This does not give the right answer at the surface
    nobs = 50
    x, y = obsgrid(-5, -5, 5, 0, nobs)
    @time Usurface, Ssurface = constdispstress(slip2dispstress, x, y, els, idx["surface"], ueff[1:2:end], ueff[2:2:end], mu, nu)
    @time Ufault, Sfault = constdispstress(slip2dispstress, x, y, els, idx["fault"], bcs[1:2:end], bcs[2:2:end], mu, nu)
    uxbem = -Usurface[:, 1] + Ufault[:, 1]
    uybem = -Usurface[:, 2] + Ufault[:, 2]

    #! Okada interior evaluation
    uxokadainterior, uyokadainterior = okadalocalinterior(x, y)

    #! 18 panel plot
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    nrows = 3
    ncols = 3
    xmat = reshape(x, nobs, nobs)
    ymat = reshape(y, nobs, nobs)
    contourvec = collect(LinRange(-0.5, 0.5, 20))

    figure(figsize=(30, 20))

    # Okada results
    subplot(nrows, ncols, 1)
    quiver(x, y, uxokadainterior, uyokadainterior)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 2)
    umat = reshape(uxokadainterior, nobs, nobs)
    contourf(xmat, ymat, umat, levels=contourvec)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(xmat, ymat, umat, levels=contourvec, colors="k", linewidths=0.5)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 3)
    umat = reshape(uyokadainterior, nobs, nobs)
    contourf(xmat, ymat, umat, levels=contourvec)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(xmat, ymat, umat, levels=contourvec, colors="k", linewidths=0.5)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)


    # IBEM results
    subplot(nrows, ncols, 4)
    quiver(x, y, uxbem, uybem)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)")
    ylabel("y (m)")
    title("ux (fault)")

    subplot(nrows, ncols, 5)
    umat = reshape(uxbem, nobs, nobs)
    contourf(xmat, ymat, umat, levels=contourvec)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(xmat, ymat, umat, levels=contourvec, colors="k", linewidths=0.5)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 6)
    umat = reshape(uybem, nobs, nobs)
    contourf(xmat, ymat, umat, levels=contourvec)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(xmat, ymat, umat, levels=contourvec, colors="k", linewidths=0.5)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    # Residuals
    subplot(nrows, ncols, 7)
    quiver(x, y, uxbem .- uxokadainterior, uybem - uyokadainterior)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 8)
    umat = reshape(uxbem .- uxokadainterior, nobs, nobs)
    contourf(xmat, ymat, umat, levels=contourvec)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(xmat, ymat, umat, levels=contourvec, colors="k", linewidths=0.5)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 9)
    umat = reshape(uybem .- uyokadainterior, nobs, nobs)
    contourf(xmat, ymat, umat, levels=contourvec)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(xmat, ymat, umat, levels=contourvec, colors="k", linewidths=0.5)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)



    @infiltrate
    return

    #! Okada solution
    xokada = collect(LinRange(-5, 5, 1000))
    uxokada, uyokada = okadalocal(xokada)

    #! Plot comparison between BEM and Okada
    figure(figsize = (20, 20))

    ax = subplot(2, 1, 1)
    plot(xokada, uxokada, "-k", linewidth=linewidth, label="Okada (halfspace)")
    plot(xbem, uhalfspacehacky[1:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label="hacky BEM (halfspace)")
    plot(xbem, uhalfspaceindirect[1:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label="indirect BEM (halfspace)")
    plot(xbem, uinterior[1:2:end], ".g", markeredgewidth=linewidth, markersize=markersize, label = "interior BEM (halfspace)")
    ylabel(L"$u_x$ (m)", fontsize=fontsize)
    plotformat(fontsize)

    ax = subplot(2, 1, 2)
    plot(xokada, uyokada, "-k", linewidth=linewidth, label="Okada (halfspace)")
    plot(xbem, uhalfspacehacky[2:2:end], "bx", markeredgewidth=linewidth, markersize=markersize, label = "hacky BEM (halfspace)")
    plot(xbem, uhalfspaceindirect[2:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "indirect BEM (halfspace)")
    plot(xbem, uinterior[2:2:end], ".g", markeredgewidth=linewidth, markersize=markersize, label = "interior BEM (halfspace)")
    ylabel(L"$u_y$ (m)", fontsize=fontsize)
    plotformat(fontsize)

    show()
    @infiltrate
end
okadaindirect()
