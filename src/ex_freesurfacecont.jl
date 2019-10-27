using Revise
using PyCall
using PyPlot
using Bem2d


function ex_freesurface()
    μ = 30e9
    ν = 0.25

    # Free surface
    els = Elements(Int(1e5))
    nfreesurf = 20
    x1, y1, x2, y2 = discretizedline(-5, 0, 5, 0, nfreesurf)

    for i in 1:length(x1)
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
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
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

    # Okada solution
    ow = pyimport("okada_wrapper")# from okada_wrapper import dc3dwrapper
    xokada = collect(LinRange(-5, 5, 10000))
    yokada = -0.1 * ones(size(xokada))
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


    # Off-fault stresses
    faultidx = findall(x->x == "fault", els.name)
    freesurfidx = findall(x->x == "freesurf", els.name)
    ufaultconstvol, σfaultconstvol = constuσ(slip2uσ, xokada, yokada, els, faultidx, faultslipconst[1:2:end], faultslipconst[2:2:end], μ, ν)
    ufreesurfaceconstvol, σfreesurfaceconstvol = constuσ(slip2uσ, xokada, yokada, els, freesurfidx, ufreesurfaceconst[1:2:end], ufreesurfaceconst[2:2:end], μ, ν)
    σconst = σfaultconstvol + σfreesurfaceconstvol

    qux = transpose(reshape(ufreesurfacequad[1:2:end], 3, nfreesurf))
    quy = transpose(reshape(ufreesurfacequad[2:2:end], 3, nfreesurf))
    ufaultquadvol, σfaultquadvol = quaduσ(slip2uσ, xokada, yokada, els, faultidx, transpose(faultslipquad[1:2:end]), transpose(faultslipquad[2:2:end]), μ, ν)
    ufreesurfacequadvol, σfreesurfacequadvol = quaduσ(slip2uσ, xokada, yokada, els, freesurfidx, qux, quy, μ, ν)
    σquad = σfaultquadvol + σfreesurfacequadvol

    # Plot ux and uy profiles
    fontsize = 24
    markersize = 15
    linewidth = 2.0
    close("all")
    figure(figsize = (15, 15))

    ax = subplot(3, 1, 1)
    plot(xokada, log10.(abs.(σconst[:, 1])), "-r", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σquad[:, 1])), "-b", linewidth=linewidth, label="3QN BEM")
    gca().set_xlim([-5, 5]);
    gca().set_ylim([6, 12])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([6, 9, 12])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{xx}$ (Pa)", fontsize=fontsize)

    ax = subplot(3, 1, 2)
    plot(xokada, log10.(abs.(σconst[:, 2])), "-r", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σquad[:, 2])), "-b", linewidth=linewidth, label="3QN BEM")
    gca().set_xlim([-5, 5]);
    gca().set_ylim([6, 12])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([6, 9, 12])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{yy}$ (Pa)", fontsize=fontsize)

    ax = subplot(3, 1, 3)
    plot(xokada, log10.(abs.(σconst[:, 3])), "-r", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σquad[:, 3])), "-b", linewidth=linewidth, label="3QN BEM")
    gca().set_xlim([-5, 5]);
    gca().set_ylim([6, 12])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([6, 9, 12])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{xy}$ (Pa)", fontsize=fontsize)

    tight_layout()
    show()

end
ex_freesurface()
