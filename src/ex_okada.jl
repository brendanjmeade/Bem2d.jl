using Revise
using PyCall
using PyPlot
using LinearAlgebra
using Bem2d

function plotfields_contours_local(els, xobs, yobs, idx, field, title)
    ncontours = 40
    xlim = [minimum(xobs) maximum(xobs)]
    ylim = [minimum(yobs) maximum(yobs)]
    subplot(2, 3, idx)
    scale = 1.0
    fieldmax = maximum(@.abs(field))
    contourf(xobs, yobs, reshape(field, size(xobs)), ncontours,
        vmin = -scale * fieldmax, vmax = scale * fieldmax, cmap = rycroftcmap())
    clim(-scale * fieldmax, scale * fieldmax)
    colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    contour(xobs, yobs, reshape(field, size(xobs)), ncontours,
        vmin = -scale * fieldmax, vmax = scale * fieldmax, linewidths = 0.25, colors = "k")
    PyPlot.title(title)
    # stylesubplots(xlim, ylim)
    plotelements(els)
    return nothing
end

function plotfields_local(els, xobs, yobs, disp, stress, suptitle)
    figure(figsize = (25, 15))
    subplot(2, 3, 1)
    quiver(xobs[:], yobs[:], disp[:, 1], disp[:, 2], units = "width", color = "b")
    # stylesubplots([minimum(xobs) maximum(xobs)], [minimum(yobs) maximum(yobs)])
    Bem2d.plotelements(els)
    PyPlot.title("displacements")

    contourvec = collect(LinRange(-0.5, 0.5, 50))
    subplot(2, 3, 2)
    contourf(xobs, yobs, reshape(disp[:, 1], size(xobs)), contourvec, cmap = rycroftcmap())
    colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    contour(xobs, yobs, reshape(disp[:, 1], size(xobs)), contourvec, linewidths = 0.25, colors = "k")
    PyPlot.title(L"u_x")
    plotelements(els)

    subplot(2, 3, 3)
    contourf(xobs, yobs, reshape(disp[:, 2], size(xobs)), contourvec, cmap = rycroftcmap())
    colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    contour(xobs, yobs, reshape(disp[:, 2], size(xobs)), contourvec, linewidths = 0.25, colors = "k")
    PyPlot.title(L"u_y")
    plotelements(els)

    plotfields_contours_local(els, xobs, yobs, 4, stress[:, 1], L"\sigma_{xx}")
    plotfields_contours_local(els, xobs, yobs, 5, stress[:, 2], L"\sigma_{yy}")
    plotfields_contours_local(els, xobs, yobs, 6, stress[:, 3], L"\sigma_{xy}")
    PyPlot.suptitle(suptitle, fontsize=20)
    show()
    return nothing
end

function ex_okada()
    mu = 30e9
    nu = 0.25

    # Flat fault
    nfault = 1
    x1, y1, x2, y2 = Bem2d.discretizedline(-0.5, 0, 0.5, 0, nfault)
    els = Bem2d.Elements(Int(1e5))
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    obswidth = 3
    npts = 99
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Constant slip fault
    UfUss, σfUss = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["fault"], 1, 0, mu, nu)
    UfUts, σfUts = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["fault"], 0, 1, mu, nu)

    # Okada solution
    ow = pyimport("okada_wrapper") # from okada_wrapper import dc3dwrapper
    uokadass = zeros(length(x), 2)
    σokadass = zeros(length(x), 3)
    uokadats = zeros(length(x), 2)
    σokadats = zeros(length(x), 3)
    deep = 10000.0
    for i in 1:length(x)
        _, u, ∇u = ow.dc3dwrapper(0.6, [0.0, x[i], y[i]-deep], 0.0 + deep, 0.0, [-1000000, 1000000], [-0.5, 0.5], [0.0, 1.0, 0.0])
        ∇u = [∇u[2,2] ∇u[2,3] ; ∇u[3,2] ∇u[3,3]]
        ϵ = @. 0.5 * (∇u + transpose(∇u))
        σ = mu*I(2)*tr(ϵ) + 2*mu*ϵ
        uokadass[i, 1] = u[2]
        uokadass[i, 2] = u[3]
        σokadass[i, 1] = σ[1, 1]
        σokadass[i, 2] = σ[2, 2]
        σokadass[i, 3] = σ[1, 2]

        _, u, ∇u = ow.dc3dwrapper(0.6, [0.0, x[i], y[i]-deep], 0.0 + deep, 0.0, [-1000000, 1000000], [-0.5, 0.5], [0.0, 0.0, 1.0])
        ∇u = [∇u[2,2] ∇u[2,3] ; ∇u[3,2] ∇u[3,3]]
        ϵ = @. 0.5 * (∇u + transpose(∇u))
        σ = mu*I(2)*tr(ϵ) + 2*mu*ϵ
        uokadats[i, 1] = u[2]
        uokadats[i, 2] = u[3]
        σokadats[i, 1] = σ[1, 1]
        σokadats[i, 2] = σ[2, 2]
        σokadats[i, 3] = σ[1, 2]
    end

    PyPlot.close("all")
    plotfields_local(els, reshape(x, npts, npts), reshape(y, npts, npts), uokadass, σokadass, "Okada (strike-slip)")
    plotfields_local(els, reshape(x, npts, npts), reshape(y, npts, npts), UfUss, σfUss, "BEM (strike-slip)")
    plotfields_local(els, reshape(x, npts, npts), reshape(y, npts, npts), uokadass - UfUss, σokadass - σfUss, "Residuals (strike-slip)")
    plotfields_local(els, reshape(x, npts, npts), reshape(y, npts, npts), uokadats, σokadats, "Okada (tensile-slip)")
    plotfields_local(els, reshape(x, npts, npts), reshape(y, npts, npts), UfUts, σfUts, "BEM (tensile-slip)")
    plotfields_local(els, reshape(x, npts, npts), reshape(y, npts, npts), uokadats - UfUts, σokadats - σfUts, "Residuals (tensile-slip)")
end
ex_okada()
