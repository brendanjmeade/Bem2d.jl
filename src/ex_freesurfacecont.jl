using Revise
using PyCall
using PyPlot
using Bem2d


function ex_freesurface()
    μ = 30e9
    ν = 0.25

    # Free surface
    els = Elements(Int(1e5))
    # nfreesurf = 200
    # x1, y1, x2, y2 = discretizedline(-5, 0, 5, 0, nfreesurf)
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-50, 0, 50, 0, nfreesurf)
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
    ∂u1const, ∂σ1const, ∂t1const = ∂constuσ(slip2uσ, els, getidx("fault", els), getidx("freesurf", els), μ, ν)
    ∂u2const, ∂σ2const, ∂t2const = ∂constuσ(slip2uσ, els, getidx("freesurf", els), getidx("freesurf", els), μ, ν)
    faultslipconst = sqrt(2) / 2 * [1 ; 1]
    ufullspaceconst = ∂u1const * faultslipconst
    ufreesurfaceconst = inv(∂t2const) * (∂t1const * faultslipconst)
    xplotconst = els.xcenter[getidx("freesurf", els)]

    # Quadratic slip fault
    ∂u1quad, ∂σ1quad, ∂t1quad = ∂quaduσ(slip2uσ, els, getidx("fault", els), getidx("freesurf", els), μ, ν)
    ∂u2quad, ∂σ2quad, ∂t2quad = ∂quaduσ(slip2uσ, els, getidx("freesurf", els), getidx("freesurf", els), μ, ν)
    xplotquad = sort(els.xnodes[getidx("freesurf", els), :][:])
    faultslipquad = sqrt(2) ./ 2 * ones(6)
    ufullspacequad = ∂u1quad * faultslipquad
    ufreesurfacequad = inv(∂t2quad) * (∂t1quad * faultslipquad)

    # Okada solution
    ow = pyimport("okada_wrapper")# from okada_wrapper import dc3dwrapper
    xokada = collect(LinRange(-5, 5, 1000))
    # yokada = -0.1 * ones(size(xokada))
    yokada = -1.8e-1 * ones(size(xokada))

    uxokada = zeros(length(xokada))
    uyokada = zeros(length(xokada))
    σxxokada = zeros(length(xokada))
    σyyokada = zeros(length(xokada))
    σxyokada = zeros(length(xokada))

    for i in 1:length(xokada)
        # Fault dipping at 45 degrees
        _, u, s = ow.dc3dwrapper(
            2.0 / 3.0,
            [0, xokada[i] + 0.5, yokada[i]],
            0.5,
            45,
            [-1e5, 1e5],
            [-sqrt(2) / 2, sqrt(2) / 2],
            [0.0, 1.0, 0.0],
        )
        uxokada[i] = u[2]
        uyokada[i] = u[3]
        dgtxx = s[2, 2]
        dgtyy = s[3, 3]
        dgtxy = s[3, 2]
        dgtyx = s[2, 3]
        exx = dgtxx
        eyy = dgtyy
        exy = 0.5 * (dgtyx + dgtxy)
        sxx = μ * (exx + eyy) + 2 * μ * exx
        syy = μ * (exx + eyy) + 2 * μ * eyy
        sxy = 2 * μ * exy
        σxxokada[i] = sxx
        σyyokada[i] = syy
        σxyokada[i] = sxy
    end

    # Off-fault stresses
    faultidx = findall(x->x == "fault", els.name)
    freesurfidx = findall(x->x == "freesurf", els.name)

    ufaultconstvol, σfaultconstvol = constuσ(slip2uσ, xokada, yokada, els, faultidx, faultslipconst[1:2:end], faultslipconst[2:2:end], μ, ν)
    ufreesurfaceconstvol, σfreesurfaceconstvol = constuσ(slip2uσ, xokada, yokada, els, freesurfidx, ufreesurfaceconst[1:2:end], ufreesurfaceconst[2:2:end], μ, ν)
    uconst = ufaultconstvol - ufreesurfaceconstvol # Note negative sign
    σconst = σfaultconstvol - σfreesurfaceconstvol # Note negative sign

    qux = transpose(reshape(ufreesurfacequad[1:2:end], 3, nfreesurf))
    quy = transpose(reshape(ufreesurfacequad[2:2:end], 3, nfreesurf))
    ufaultquadvol, σfaultquadvol = quaduσ(slip2uσ, xokada, yokada, els, faultidx, transpose(faultslipquad[1:2:end]), transpose(faultslipquad[2:2:end]), μ, ν)
    ufreesurfacequadvol, σfreesurfacequadvol = quaduσ(slip2uσ, xokada, yokada, els, freesurfidx, qux, quy, μ, ν)
    uquad = ufaultquadvol - ufreesurfacequadvol # Note negative sign
    σquad = σfaultquadvol - σfreesurfacequadvol # Note negative sign

    # Plot ux and uy profiles
    fontsize = 24
    markersize = 15
    linewidth = 5.0
    close("all")
    figure(figsize = (30, 15))

    # ax = subplot(2, 3, 2)
    # plot(xokada, uconst[:, 1], "-c", linewidth=linewidth, label = "CS BEM")
    # plot(xokada, uquad[:, 1], "-r", linewidth=linewidth, label = "3QN BEM")
    # plot(xokada, uxokada, "--b", linewidth=2.0, label="Okada")
    # gca().set_xlim([-5, 5])
    # gca().set_ylim([-1.0, 1.0])
    # gca().set_xticks([-5, 0, 5])
    # gca().set_yticks([-1, 0, 1])
    # legend(fontsize=fontsize)
    # ax.tick_params("both", labelsize = fontsize)
    # xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$u_x$ (m)", fontsize=fontsize)
    #
    # ax = subplot(2, 3, 3)
    # plot(xokada, uconst[:, 2], "-c", linewidth=linewidth, label = "CS BEM")
    # plot(xokada, uquad[:, 2], "-r", linewidth=linewidth, label = "3QN BEM")
    # plot(xokada, uyokada, "--b", linewidth=2.0, label="Okada")
    # gca().set_xlim([-5, 5])
    # gca().set_ylim([-1.0, 1.0])
    # # gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    # gca().set_xticks([-5, 0, 5])
    # gca().set_yticks([-1, 0, 1])
    # legend(fontsize=fontsize, frameon=false)
    # ax.tick_params("both", labelsize = fontsize)
    # xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$u_y$ (m)", fontsize=fontsize)

    ax = subplot(2, 3, 1)
    plot(xokada, log10.(abs.(σconst[:, 1])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σquad[:, 1])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(σxxokada[:, 1])), "--k", linewidth=2.0, label="Okada")
    gca().set_xlim([-5, 5])
    gca().set_ylim([6, 12])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([6, 9, 12])
    legend(fontsize=fontsize, frameon=false)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{xx}$ (Pa)", fontsize=fontsize)

    ax = subplot(2, 3, 2)
    plot(xokada, log10.(abs.(σconst[:, 2])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σquad[:, 2])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(σyyokada[:, 1])), "--k", linewidth=2.0, label="Okada")
    gca().set_xlim([-5, 5])
    gca().set_ylim([6, 12])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([6, 9, 12])
    legend(fontsize=fontsize, frameon=false)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{yy}$ (Pa)", fontsize=fontsize)

    ax = subplot(2, 3, 3)
    plot(xokada, log10.(abs.(σconst[:, 3])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σquad[:, 3])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(σxyokada[:, 1])), "--k", linewidth=2.0, label="Okada")
    gca().set_xlim([-5, 5])
    gca().set_ylim([6, 12])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([6, 9, 12])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{xy}$ (Pa)", fontsize=fontsize)

    ax = subplot(2, 3, 4)
    σxxconsterror = 100 * (σconst[:, 1] - σxxokada[:, 1]) ./ σxxokada[:, 1]
    σxxquaderror = 100 * (σquad[:, 1] - σxxokada[:, 1]) ./ σxxokada[:, 1]
    plot(xokada, log10.(abs.(σxxconsterror)), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σxxquaderror)), "-r", linewidth=linewidth, label="3QN BEM")
    gca().set_xlim([-5, 5])
    gca().set_ylim([-3, 6])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsize, frameon=false)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{xx} \, \%$ error", fontsize=fontsize)

    ax = subplot(2, 3, 5)
    σyyconsterror = 100 * (σconst[:, 2] - σyyokada[:, 1]) ./ σyyokada[:, 1]
    σyyquaderror = 100 * (σquad[:, 2] - σyyokada[:, 1]) ./ σyyokada[:, 1]
    plot(xokada, log10.(abs.(σyyconsterror)), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σyyquaderror)), "-r", linewidth=linewidth, label="3QN BEM")
    gca().set_xlim([-5, 5])
    gca().set_ylim([-3, 6])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsize, frameon=false)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{yy} \, \%$ error", fontsize=fontsize)

    ax = subplot(2, 3, 6)
    σxyconsterror = 100 * (σconst[:, 3] - σxyokada[:, 1]) ./ σxyokada[:, 1]
    σxyquaderror = 100 * (σquad[:, 3] - σxyokada[:, 1]) ./ σxyokada[:, 1]
    plot(xokada, log10.(abs.(σxyconsterror)), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(σxyquaderror)), "-r", linewidth=linewidth, label="3QN BEM")
    gca().set_xlim([-5, 5])
    gca().set_ylim([-3, 6])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsize, frameon=false)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{xy} \, \%$ error", fontsize=fontsize)

    tight_layout()
    show()

end
ex_freesurface()
