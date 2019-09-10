using Revise
using PyCall
using PyPlot
using Bem2d

function ex_constquadpartials()
    # Material and geometric constants
    μ = 3e10
    ν = 0.25
    nels = 50
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
    close("all")
    # matshow(∂uconst); title("∂uconst"); colorbar()
    # matshow(∂σconst); title("∂σconst"); colorbar()
    # matshow(∂tconst); title("∂tconst"); colorbar()
    # matshow(∂uquad); title("∂uquad"); colorbar()
    # matshow(∂σquad); title("∂σquad"); colorbar()
    # matshow(∂tquad); title("∂tquad"); colorbar()

    # Evaluation points and slip
    xcenters = unique(els.xcenter[1:els.endidx])
    ycenters = unique(els.ycenter[1:els.endidx])
    xnodes = unique(els.xnodes[1:els.endidx, :])
    ynodes = unique(els.ynodes[1:els.endidx, :])
    slipconst = zeros(2 * nels)
    slipquad = zeros(6 * nels)

    # Strike-slip
    slipquad[1:2:end] .= 1  # constant x-slip global
    slipconst[1:2:end] .= 1  # constant x-slip global
    # Tensile-slip
    # slipquad[2:2:end] .= 1  # constant x-slip global
    # slipconst[2:2:end] .= 1  # constant x-slip global

    # Predict on fault displacements, stresses, and tractions
    uconst = ∂uconst * slipconst
    σconst = ∂σconst * slipconst
    tconst = ∂tconst * slipconst
    uquad = ∂uquad * slipquad
    σquad = ∂σquad * slipquad
    tquad = ∂tquad * slipquad

    # Plot geometry of elements
    figure(figsize = (15, 10))
    subplot(2, 3, 1)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-b", color = "b", linewidth = 0.5)
        plot([els.x1[i] els.x2[i]], [els.y1[i] els.y2[i]], ".r", markersize = 10, linewidth = 0.5)
        text(els.xcenter[i], els.ycenter[i], string(i), horizontalalignment = "center", verticalalignment = "center", fontsize = 8)
    end
    ylabel("y (m)"); xlim([-10000, 10000]); xticks([-10000, 0, 10000])

    subplot(2, 3, 2)
    plot(xcenters, uconst[1:2:end], "ob", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$u_x$ constant")
    plot(xnodes, uquad[1:2:end], "+r", markeredgewidth = 0.5, label = L"$u_x$ quadratic")
    legend(loc = "lower right"); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ylabel("displacement (m)")

    subplot(2, 3, 3)
    plot(xcenters, uconst[2:2:end], "ob", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$u_y$ constant")
    plot(xnodes, uquad[2:2:end], "+r", markeredgewidth = 0.5, label = L"$u_y$ quadratic")
    legend(loc = "lower right"); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    ylabel("displacement (m)")

    subplot(2, 3, 4)
    plot(xcenters, σconst[1:3:end], "ob", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{xy}$ constant")
    plot(xnodes, σquad[1:3:end], "+r", markeredgewidth = 0.5, label = L"$σ_{xx}$ quadratic")
    legend(loc = "lower right"); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel("x (m)"); ylabel("stress (Pa)")

    subplot(2, 3, 5)
    plot(xcenters, σconst[2:3:end], "ob", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{yy}$ constant")
    plot(xnodes, σquad[2:3:end], "+r", markeredgewidth = 0.5, label = L"$σ_{yy}$ quadratic")
    legend(loc = "lower right"); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel("x (m)"); ylabel("stress (Pa)")

    subplot(2, 3, 6)
    plot(xcenters, σconst[3:3:end], "ob", markerfacecolor = "none", markeredgewidth = 0.5, label = L"$σ_{xy}$ constant")
    plot(xnodes, σquad[3:3:end], "+r", markeredgewidth = 0.5, label = L"$σ_{xy}$ quadratic")
    legend(loc = "lower right"); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel("x (m)"); ylabel("stress (Pa)")
    show()
end
ex_constquadpartials()
