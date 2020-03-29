using Revise
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function stylesubplots_local()
    gca().set_aspect("equal")
    xlim([-200, 200])
    ylim([-200, 200])
    return nothing
end

function subplotlocal(nrows, ncols, plotidx, els, x, y, field, contours, titlestring)
    fontsize = 30
    cmap = get_cmap("seismic")
    subplot(nrows, ncols, plotidx)
    contourf(x, y, reshape(field, size(x)), contours, cmap=cmap)
    colorbar(fraction=0.020, pad=0.05, extend="both")
    contour(x, y, reshape(field, size(x)), contours, linewidths=0.25, colors="k")
    stylesubplots_local()
    plotelements(els)
    title(titlestring, fontsize=fontsize)
    return nothing
end

function ex_flamantfinite()
    obswidth = 200

    # Uniform mesh spacing
    npts = 51
    x1, y1, x2, y2 = discretizedline(-2000, 0, 2000, 0, npts)
    @show a = x1[2] - x1[1]
    mididx = Int(floor(length(x1)/2)+1)
    
    mu = 3e10
    nu = 0.25
    npts = 50
    x, y = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    nels = length(x1)
    P = 1.0
    fx = zeros(nels)
    fy = zeros(nels)
    fy[mididx] = 1.0
    
    #! Rescale forcing for smaller central element. This is strange
    fxscaled = @. fx / (2*a) #! Strange length dependent normalization
    fyscaled = @. fy / (2*a) #! Strange length dependent normalization

    #! Analytic finite Flamant solution
    Danalytic = zeros(length(x), 2)
    Sanalytic = zeros(length(x), 3)
    r1 = @. (x-a)^2 + y^2
    r2 = @. (x+a)^2 + y^2
    theta1 = @. atan(y, x-a)
    theta2 = @. atan(y, x+a)
    Sanalytic[:, 1] = @. -P/pi * (theta1-theta2 + y*(x-a)*r1^-2 - y*(x+a)*r2^-2)
    Sanalytic[:, 2] = @. -P/pi * (theta1-theta2 - y*(x-a)*r1^-2 + y*(x+a)*r2^-2)
    Sanalytic[:, 3] = @. -P/pi * y^2 * (r1^-2 - r2^-2)

    #! BEM solution
    els = Elements(Int(1e5))
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "line"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Kernels
    Tstar, Sstar, Hstar = partialsconstdispstress(slip2dispstress, els, idx["line"], idx["line"], mu, nu)
    Ustar, Dstar, Astar = partialsconstdispstress(trac2dispstress, els, idx["line"], idx["line"], mu, nu)
    
    #! Direct BEM
    dispall = (inv(Tstar + 0.5 * I(size(Tstar)[1]))) * Ustar * interleave(fxscaled, fyscaled)

    #! Try Indirect DDM...because it worked for the Brazil test
    # dispall = -Ustar * interleave(fx, fy)

    #! Streses from tractions and induced displacements
    Dtrac, Strac = constdispstress(trac2dispstress, x, y, els, idx["line"], fxscaled, fyscaled, mu, nu)
    Ddisp, Sdisp = constdispstress(slip2dispstress, x, y, els, idx["line"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    #! Set everything in the upper half-plane to zero because that is the exterior solution
    to_nan_idx = findall(x -> x > 0.0, y)
    Danalytic[to_nan_idx, :] .= NaN
    Sanalytic[to_nan_idx, :] .= NaN
    Dtrac[to_nan_idx, :] .= NaN
    Ddisp[to_nan_idx, :] .= NaN
    Strac[to_nan_idx, :] .= NaN
    Sdisp[to_nan_idx, :] .= NaN

    #! Multi-panel plot
    close("all")
    # Dbem = @. Dtrac - Ddisp # Direct BEM
    # Sbem = @. Strac - Sdisp # Direct BEM
    Dbem = Ddisp # DDM
    Sbem = Sdisp # DDM
    xplot = reshape(x, npts, npts)
    yplot = reshape(y, npts, npts)
    contoursdisp = 51
    contourvecstress = 51
    cmap = get_cmap("seismic")
    fontsize = 30
    nrows = 3
    ncols = 6
    figure(figsize=(30, 20))

    #! Analytic solution
    subplot(nrows, ncols, 1)
    quiver(x[:], y[:], Danalytic[:, 1], Danalytic[:, 2], units="width", color="b")
    stylesubplots_local()
    plotelements(els)
    title("u", fontsize=fontsize)
    subplotlocal(nrows, ncols, 2, els, xplot, yplot, Danalytic[:, 1], contoursdisp, "ux (analytic)")
    subplotlocal(nrows, ncols, 3, els, xplot, yplot, Danalytic[:, 2], contoursdisp, "uy (analytic)")
    subplotlocal(nrows, ncols, 4, els, xplot, yplot, Sanalytic[:, 1], contoursdisp, "sxx (analytic)")
    subplotlocal(nrows, ncols, 5, els, xplot, yplot, Sanalytic[:, 2], contoursdisp, "syy (analytic)")
    subplotlocal(nrows, ncols, 6, els, xplot, yplot, Sanalytic[:, 3], contoursdisp, "sxy (analytic)")

    # #! Fields from second model
    subplot(nrows, ncols, 7)
    quiver(x[:], y[:], Dbem[:, 1], Dbem[:, 2], units="width", color="b")
    stylesubplots_local()
    plotelements(els)
    title("u", fontsize=fontsize)
    subplotlocal(nrows, ncols, 8, els, xplot, yplot, Dbem[:, 1], contoursdisp, "ux (BEM)")
    subplotlocal(nrows, ncols, 9, els, xplot, yplot, Dbem[:, 2], contoursdisp, "uy (BEM)")
    subplotlocal(nrows, ncols, 10, els, xplot, yplot, Sbem[:, 1], contoursdisp, "sxx (BEM)")
    subplotlocal(nrows, ncols, 11, els, xplot, yplot, Sbem[:, 2], contoursdisp, "syy (BEM)")
    subplotlocal(nrows, ncols, 12, els, xplot, yplot, Sbem[:, 3], contoursdisp, "sxy (BEM)")

    # #! Residuals
    subplot(3, 6, 13)
    quiver(x[:], y[:], Danalytic[:, 1]-Dbem[:, 1], Danalytic[:, 2]-Dbem[:, 2], units="width", color="b")
    stylesubplots_local()
    plotelements(els)
    title("u", fontsize=fontsize)
    subplotlocal(nrows, ncols, 14, els, xplot, yplot, Danalytic[:, 1]-Dbem[:, 1], contoursdisp, "ux residuals")
    subplotlocal(nrows, ncols, 15, els, xplot, yplot, Danalytic[:, 2]-Dbem[:, 2], contoursdisp, "uy residuals")
    subplotlocal(nrows, ncols, 16, els, xplot, yplot, Sanalytic[:, 1]-Sbem[:, 1], contoursdisp, "sxx residuals")
    subplotlocal(nrows, ncols, 17, els, xplot, yplot, Sanalytic[:, 2]-Sbem[:, 2], contoursdisp, "syy residuals")
    subplotlocal(nrows, ncols, 18, els, xplot, yplot, Sanalytic[:, 3]-Sbem[:, 3], contoursdisp, "sxy residuals")

    tight_layout()
    show()
end
ex_flamantfinite()
