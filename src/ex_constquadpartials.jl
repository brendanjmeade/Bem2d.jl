using Revise
using PyCall
using PyPlot
using Bem2d

function plotpartials(els, xcenters, ycenters, xnodes, ynodes, dispconst, stressconst, tracconst, dispquad, stressquad, tracquad)
    fontsize = 6
    figure(figsize = (12, 8))
    ax = subplot(2, 3, 1)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-b", color = "b", linewidth = 0.5)
        plot([els.x1[i] els.x2[i]], [els.y1[i] els.y2[i]], ".r", markersize = 10, linewidth = 0.5)
    end
    axis("off")
    ylabel("y (m)", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ax.tick_params("both", labelsize = fontsize)

    ax = subplot(2, 3, 2)
    plot(xcenters, dispconst[1:2:end], ".b",  markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$u_x$ constant", zorder = 2)
    plot(xnodes, dispquad[1:2:end], "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$u_x$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ylabel(L"$u$ (m)", fontsize = fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylim([-0.6, 0.6])

    ax = subplot(2, 3, 3)
    plot(xcenters, dispconst[2:2:end], ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$u_y$ constant", zorder = 2)
    plot(xnodes, dispquad[2:2:end], "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$u_y$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ylabel(L"$u$ (m)", fontsize = fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylim([-0.6, 0.6])

    ax = subplot(2, 3, 4)
    plot(xcenters, stressconst[1:3:end] ./ 1e3, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{xx}$ constant", zorder = 2)
    plot(xnodes, stressquad[1:3:end] ./ 1e3, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{xx}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylim([-1e7 / 1e3, 1e7 / 1e3])

    ax = subplot(2, 3, 5)
    plot(xcenters, stressconst[2:3:end] ./ 1e3, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{yy}$ constant", zorder = 2)
    plot(xnodes, stressquad[2:3:end] ./ 1e3, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{yy}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params.("both", labelsize = fontsize)
    ylim([-1e7 / 1e3, 1e7 / 1e3])

    ax = subplot(2, 3, 6)
    plot(xcenters, stressconst[3:3:end] ./ 1e3, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{xy}$ constant", zorder = 2)
    plot(xnodes, stressquad[3:3:end] ./ 1e3, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{xy}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylim([-1e7 / 1e3, 1e7 / 1e3])
    tight_layout(); show()
end

function ex_constquadpartials()
    # Material and geometric constants
    mu = 3e10
    nu = 0.25
    nels = 5
    els = Elements(Int(1e5))
    L = 5000
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Convenience structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    partialsquad = initpartials(els)

    # Partial derivatves
    partialsconst["disp"]["fault"]["fault"], partialsconst["stress"]["fault"]["fault"], partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    partialsquad["disp"]["fault"]["fault"], partialsquad["stress"]["fault"]["fault"], partialsquad["trac"]["fault"]["fault"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)

    # Evaluation points and slip
    xcenters = els.xcenter[1:els.endidx]
    ycenters = els.ycenter[1:els.endidx]
    xnodes = sort(els.xnodes[1:els.endidx, :][:])
    ynodes = sort(els.ynodes[1:els.endidx, :][:])

    # Parameters for the test cases
    slipconst = zeros(2, 2 * nels)
    slipquad = zeros(2, 6 * nels)

    # Constant strike-slip
    slipconst[1, 1:2:end] .= 1
    slipquad[1, 1:2:end] .= 1

    # Linear strike-slip
    slope = 0.0001
    slipconst[2, 1:2:end] = slope * xcenters
    slipquad[2, 1:2:end] = slope * xnodes

    # Predict on-fault displacements, stresses, and tractions
    close("all")
    for i in 1:size(slipconst)[1]
        dispconst = partialsconst["disp"]["fault"]["fault"] * slipconst[i, :]
        stressconst = partialsconst["stress"]["fault"]["fault"] * slipconst[i, :]
        tracconst = partialsconst["trac"]["fault"]["fault"] * slipconst[i, :]
        dispquad = partialsquad["disp"]["fault"]["fault"] * slipquad[i, :]
        stressquad = partialsquad["stress"]["fault"]["fault"] * slipquad[i, :]
        tracquad = partialsquad["trac"]["fault"]["fault"] * slipquad[i, :]
        plotpartials(els, xcenters, ycenters, xnodes, ynodes, dispconst, stressconst, tracconst, dispquad, stressquad, tracquad)
    end
end
ex_constquadpartials()
