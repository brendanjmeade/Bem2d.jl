using Revise
using PyCall
using PyPlot
using Bem2d

function ex_freesurface()
    μ = 30e9
    ν = 0.25

    # Free surface
    els = Elements(Int(1e5))
    nfreesurf = 40
    x1, y1, x2, y2 = discretizedline(-5, 0, 5, 0, nfreesurf)
    els.x1[els.endidx + 1:els.endidx + nfreesurf] = x1
    els.y1[els.endidx + 1:els.endidx + nfreesurf] = y1
    els.x2[els.endidx + 1:els.endidx + nfreesurf] = x2
    els.y2[els.endidx + 1:els.endidx + nfreesurf] = y2
    els.name[els.endidx + 1:els.endidx + nfreesurf] .= "freesurf"
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

    # Constant slip fault
    srcidx = findall(x->x == "fault", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u1const, ∂σ1const, ∂t1const = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)
    srcidx = findall(x->x == "freesurf", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u2const, ∂σ2const, ∂t2const = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)
    faultslipconst = sqrt(2) / 2 * [1 ; 1]
    ufullspaceconst = ∂u1const * faultslipconst
    ufreesurfaceconst = inv(∂t2const) * (∂t1const * faultslipconst)
    xplotconst = els.xcenter[findall(x->x == "freesurf", els.name)]

    # Quadratic slip fault
    srcidx = findall(x->x == "fault", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u1quad, ∂σ1quad, ∂t1quad = ∂quaduσ(slip2uσ, els, srcidx, obsidx, μ, ν)
    srcidx = findall(x->x == "freesurf", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u2quad, ∂σ2quad, ∂t2quad = ∂quaduσ(slip2uσ, els, srcidx, obsidx, μ, ν)
    xplotquad = sort(els.xnodes[obsidx, :][:])
    faultslipquad = sqrt(2) ./ 2 * ones(6)
    ufullspacequad = ∂u1quad * faultslipquad
    ufreesurfacequad = inv(∂t2quad) * (∂t1quad * faultslipquad)

    fontsize = 6
    markersize = 4
    linewidth = 0.5
    close("all")
    figure(figsize = (5, 5))

    ax = subplot(2, 1, 1)
    plot(xplotconst, ufullspaceconst[1:2:end], "bo", markeredgewidth=linewidth, markersize=markersize, label = "const fullspace")
    plot(xplotconst, ufreesurfaceconst[1:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    plot(xplotquad, ufullspacequad[1:2:end], "r.", markeredgewidth=linewidth, markersize=markersize, label = "quad fullspace")
    plot(xplotquad, ufreesurfacequad[1:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "quad halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-0.6, 0.6])
    gca().set_xticks([-5, 0, 5]); gca().set_yticks([-0.5, 0.0, 0.5])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$u_x$ (m)", fontsize=fontsize)

    ax = subplot(2, 1, 2)
    plot(xplotconst, ufullspaceconst[2:2:end], "bo", markeredgewidth=linewidth, markersize=markersize, label = "const fullspace")
    plot(xplotconst, ufreesurfaceconst[2:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    plot(xplotquad, ufullspacequad[2:2:end], "r.", markeredgewidth=linewidth, markersize=markersize, label = "quad fullspace")
    plot(xplotquad, ufreesurfacequad[2:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "quad halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-0.6, 0.6])
    gca().set_xticks([-5, 0, 5]); gca().set_yticks([-0.5, 0.0, 0.5])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$u_y$ (m)", fontsize=fontsize)
    show()
end
ex_freesurface()
