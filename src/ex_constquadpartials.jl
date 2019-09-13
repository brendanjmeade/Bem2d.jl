using Revise
using PyCall
using PyPlot
using Bem2d

function plotpartials(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad)
    fontsize = 6
    figure(figsize = (12, 8))
    ax = subplot(2, 3, 1)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-b", color = "b", linewidth = 0.5)
        plot([els.x1[i] els.x2[i]], [els.y1[i] els.y2[i]], ".r", markersize = 10, linewidth = 0.5)
        text(els.xcenter[i], els.ycenter[i], string(i), horizontalalignment = "center", verticalalignment = "center", fontsize = fontsize)
    end
    ylabel("y (m)", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ax[:tick_params]("both", labelsize = fontsize)

    ax = subplot(2, 3, 2)
    plot(xcenters, uconst[1:2:end], "xb", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$u_x$, constant")
    plot(xnodes, uquad[1:2:end], "sr", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$u_x$, quadratic")
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ylabel("displacement (m)", fontsize = fontsize)
    ax[:tick_params]("both", labelsize = fontsize)
    ylim([-0.6, 0.6])

    ax = subplot(2, 3, 3)
    plot(xcenters, uconst[2:2:end], "xb", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$u_y$, constant")
    plot(xnodes, uquad[2:2:end], "sr", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$u_y$, quadratic")
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ylabel("displacement (m)", fontsize = fontsize)
    ax[:tick_params]("both", labelsize = fontsize)
    ylim([-0.6, 0.6])

    ax = subplot(2, 3, 4)
    plot(xcenters, σconst[1:3:end], "xb", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{xx}$, constant")
    plot(xnodes, σquad[1:3:end], "sr", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{xx}$, quadratic")
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel("x (m)", fontsize = fontsize); ylabel("stress (Pa)", fontsize = fontsize)
    ax[:tick_params]("both", labelsize = fontsize)
    ylim([-1e7, 1e7])

    ax = subplot(2, 3, 5)
    plot(xcenters, σconst[2:3:end], "xb", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{yy}$, constant")
    plot(xnodes, σquad[2:3:end], "sr", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{yy}$, quadratic")
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel("x (m)", fontsize = fontsize); ylabel("stress (Pa)", fontsize = fontsize)
    ax[:tick_params]("both", labelsize = fontsize)
    ylim([-1e7, 1e7])

    ax = subplot(2, 3, 6)
    plot(xcenters, σconst[3:3:end], "xb", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{xy}$, constant")
    plot(xnodes, σquad[3:3:end], "sr", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{xy}$, quadratic")
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel("x (m)", fontsize = fontsize); ylabel("stress (Pa)", fontsize = fontsize)
    ax[:tick_params]("both", labelsize = fontsize)
    ylim([-1e7, 1e7])
    tight_layout(); show()
end

function ex_constquadpartials()
    # Material and geometric constants
    μ = 3e10
    ν = 0.25
    nels = 10
    els = Elements(Int(1e5))
    L = 10000
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nels)
    els.x1[els.endidx + 1:els.endidx + nels] = x1
    els.y1[els.endidx + 1:els.endidx + nels] = y1
    els.x2[els.endidx + 1:els.endidx + nels] = x2
    els.y2[els.endidx + 1:els.endidx + nels] = y2
    els.name[els.endidx + 1:els.endidx + nels] .= "fault"
    standardize_elements!(els)

    srcidx = findall(x->x == "fault", els.name)
    obsidx = findall(x->x == "fault", els.name)
    ∂uconst, ∂σconst, ∂tconst = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)
    ∂uquad, ∂σquad, ∂tquad = ∂quaduσ(slip2uσ, els, srcidx, obsidx, μ, ν)

    # Evaluation points and slip
    xcenters = els.xcenter[1:els.endidx]
    ycenters = els.ycenter[1:els.endidx]
    xnodes = sort(els.xnodes[1:els.endidx, :][:])
    ynodes = sort(els.ynodes[1:els.endidx, :][:])
    slipconst = zeros(2 * nels)
    slipquad = zeros(6 * nels)

    # Uniform strike-slip
    # slipconst[1:2:end] .= 1
    # slipquad[1:2:end] .= 1
    
    # Linear strike-slip
    slope = 0.0001
    slipconst[1:2:end] = slope * xcenters
    slipquad[1:2:end] = slope * xnodes

    # Tensile-slip
    # slipquad[2:2:end] .= 1
    # slipconst[2:2:end] .= 1

    # Predict on-fault displacements, stresses, and tractions
    uconst = ∂uconst * slipconst
    σconst = ∂σconst * slipconst
    tconst = ∂tconst * slipconst
    uquad = ∂uquad * slipquad
    σquad = ∂σquad * slipquad
    tquad = ∂tquad * slipquad

    close("all")
    plotpartials(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad)

end
ex_constquadpartials()
