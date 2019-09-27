using Revise
using PyCall
using PyPlot
using Bem2d

# From Segall equations 3.73 on page 67
function analyticthrust()
    x = collect(LinRange(-5, 5, 1000))
    δ = deg2rad(135.0)
    d = 1.0
    ξ = d / tan(δ)
    ζ = @. (x - ξ) / d
    s = 1.0
    ux = @. -s / π * (cos(δ) * (atan(ζ) - π / 2 * sign(x)) + (sin(δ) - ζ * cos(δ)) / (1 + ζ^2)^2)
    uy =  @. s / π * (sin(δ) * (atan(ζ) - π / 2 * sign(x)) + (cos(δ) + ζ * sin(δ)) / (1 + ζ^2)^2)
    return x, ux, uy
end

function ex_freesurface()
    μ = 30e9
    ν = 0.25
    println(rand(1))

    # Free surface
    els = Elements(Int(1e5))
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-50, 0, 50, 0, nfreesurf)
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

    # Analytic solution
    xanalytic, uxanalytic, uyanalytic = analyticthrust()

    # Okada solution
    ow = pyimport("okada_wrapper")# from okada_wrapper import dc3dwrapper
    xokada = collect(LinRange(-5, 5, 1000))
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

    fontsize = 6
    markersize = 4
    linewidth = 0.5
    close("all")
    figure(figsize = (3, 4))

    ax = subplot(2, 1, 1)
    # plot(xanalytic, uxanalytic, "-k", linewidth=linewidth, label="analytic")
    plot(xokada, uxokada, "-k", linewidth=linewidth, label="Okada")

    # plot(xplotconst, ufullspaceconst[1:2:end], "bo", markeredgewidth=linewidth, markersize=markersize, label = "const fullspace")
    plot(xplotconst, ufreesurfaceconst[1:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    # plot(xplotquad, ufullspacequad[1:2:end], "r.", markeredgewidth=linewidth, markersize=markersize, label = "quad fullspace")
    plot(xplotquad, ufreesurfacequad[1:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "quad halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([-1.00, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$u_x$ (m)", fontsize=fontsize)

    ax = subplot(2, 1, 2)
    # plot(xanalytic, uyanalytic, "-k", linewidth=linewidth, label="analytic")
    plot(xokada, uyokada, "-k", linewidth=linewidth, label="Okada")

    # plot(xplotconst, ufreesurfaceconst[1:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")

    # plot(xplotconst, ufullspaceconst[2:2:end], "bo", markeredgewidth=linewidth, markersize=markersize, label = "const fullspace")
    plot(xplotconst, ufreesurfaceconst[2:2:end], "b+", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    # plot(xplotquad, ufullspacequad[2:2:end], "r.", markeredgewidth=linewidth, markersize=markersize, label = "quad fullspace")
    plot(xplotquad, ufreesurfacequad[2:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "quad halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([-1.00, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$u_y$ (m)", fontsize=fontsize)
    show()
end
ex_freesurface()
