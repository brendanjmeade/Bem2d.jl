using Revise
using PyCall
using PyPlot
using Bem2d

function fig_halfspace()
    mu = 30e9
    nu = 0.25

    # Free surface
    els = Elements(Int(1e5))
    nfreesurf = 120
    x1, y1, x2, y2 = discretizedline(-50, 0, 50, 0, nfreesurf)
    for i in length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "freesurf"
    end
    standardize_elements!(els)

    # 45 degree dipping fault
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-1, -1, 0, 0, nfault)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Convenience data structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    partialsquad = initpartials(els)

    # Constant slip fault
    partialsconst["disp"]["fault"]["freesurf"], _, partialsconst["trac"]["fault"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    partialsconst["disp"]["freesurf"]["freesurf"], _, partialsconst["trac"]["freesurf"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    faultslipconst = sqrt(2) / 2 * [1 ; 1]
    dispfullspaceconst = partialsconst["disp"]["fault"]["freesurf"] * faultslipconst
    dispfreesurfaceconst = inv(partialsconst["trac"]["freesurf"]["freesurf"]) * (partialsconst["trac"]["fault"]["freesurf"] * faultslipconst)
    xplotconst = els.xcenter[idx["freesurf"]]

    # Quadratic slip fault
    partialsquad["disp"]["fault"]["freesurf"], _, partialsquad["trac"]["fault"]["freesurf"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    partialsquad["disp"]["freesurf"]["freesurf"], _, partialsquad["trac"]["freesurf"]["freesurf"] = partialsquaddispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    xplotquad = sort(els.xnodes[idx["freesurf"], :][:])
    faultslipquad = sqrt(2) ./ 2 * ones(6)
    dispfullspacequad = partialsquad["disp"]["fault"]["freesurf"] * faultslipquad
    dispfreesurfacequad = inv(partialsquad["trac"]["freesurf"]["freesurf"]) * (partialsquad["trac"]["fault"]["freesurf"] * faultslipquad)

    # Okada solution
    okadawrapper = pyimport("okada_wrapper")# from okada_wrapper import dc3dwrapper
    xokada = collect(LinRange(-5, 5, 1000))
    dispxokada = zeros(length(xokada))
    dispyokada = zeros(length(xokada))
    for i in 1:length(xokada)
        _, u, _ = okadawrapper.dc3dwrapper(2.0 / 3.0, [0, xokada[i] + 0.5, 0],
            0.5, 45, [-10000, 10000], [-sqrt(2) / 2, sqrt(2) / 2], [0.0, 1.0, 0.0])
        dispxokada[i] = u[2]
        dispyokada[i] = u[3]
    end

    fontsize = 6
    markersize = 4
    linewidth = 0.5
    close("all")
    figure(figsize = (5, 7))

    ax = subplot(2, 1, 1)
    plot(xokada, dispxokada, "-r", linewidth=linewidth, label="Okada", zorder=1)
    plot(xplotconst, dispfreesurfaceconst[1:2:end], ".k", markeredgewidth=linewidth, markersize=markersize, label = "CS BEM", zorder=3)
    plot(xplotquad, dispfreesurfacequad[1:2:end], "ok", markerfacecolor = "lightgray", markeredgewidth = 0.25, markersize=markersize, label = "3NQ BEM", zorder=2)
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, 0, 5]); gca().set_yticks([-1.00, 0.00, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$u_x$ (m)", fontsize=fontsize)

    ax = subplot(2, 1, 2)
    plot(xokada, dispyokada, "-r", linewidth=linewidth, label="Okada", zorder = 1)
    plot(xplotconst, dispfreesurfaceconst[2:2:end], ".k", markeredgewidth=linewidth, markersize=markersize, label = "CS BEM", zorder = 3)
    plot(xplotquad, dispfreesurfacequad[2:2:end], "ok", markerfacecolor = "lightgray", markeredgewidth = 0.25, markersize=markersize, label = "3NQ BEM", zorder = 2)
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, 0, 5]); gca().set_yticks([-1.00, 0.00, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$u_y$ (m)", fontsize=fontsize)
    tight_layout()
    show()
end
fig_halfspace()
