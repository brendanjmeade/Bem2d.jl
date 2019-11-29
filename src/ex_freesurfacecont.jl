using Revise
using Printf
using PyCall
using PyPlot
using Statistics
using Bem2d

function ex_freesurface()
    mu = 30e9
    nu = 0.25

    # Free surface
    els = Elements(Int(1e5))
    # nfreesurf = 20
    # x1, y1, x2, y2 = discretizedline(-50e3, 0, 50e3, 0, nfreesurf)
    nfreesurf = 100
    x1, y1, x2, y2 = discretizedline(-100e3, 0, 100e3, 0, nfreesurf)
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
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Convience data structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    partialsquad = initpartials(els)
    # _, _, partialsconst["trac"]["fault"]["freesurfflat"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurfflat"], mu, nu)

    # Constant slip fault
    ∂u1const, _, ∂t1const = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    ∂u2const, _, ∂t2const = partialsconstdispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    faultslipconst = sqrt(2) / 2 * [1 ; 1]
    dispfullspaceconst = ∂u1const * faultslipconst
    dispfreesurfaceconst = inv(∂t2const) * (∂t1const * faultslipconst)
    xplotconst = els.xcenter[idx["freesurf"]]

    # Quadratic slip fault
    ∂u1quad, _, ∂t1quad = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    ∂u2quad, _, ∂t2quad = partialsquaddispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    xplotquad = sort(els.xnodes[idx["freesurf"], :][:])
    faultslipquad = sqrt(2) ./ 2 * ones(6)
    dispfullspacequad = ∂u1quad * faultslipquad
    dispfreesurfacequad = inv(∂t2quad) * (∂t1quad * faultslipquad)

    # Okada solution
    ow = pyimport("okada_wrapper")# from okada_wrapper import dc3dwrapper
    xokada = collect(LinRange(-50e3, 50e3, 1000))
    offset = abs(els.xnodes[1, 2] - els.xnodes[1, 1])
    yokada = -offset * ones(size(xokada))
    dispxokada = zeros(length(xokada))
    dispyokada = zeros(length(xokada))
    stressxxokada = zeros(length(xokada))
    stressyyokada = zeros(length(xokada))
    stressxyokada = zeros(length(xokada))

    for i in 1:length(xokada)
        # Fault dipping at 45 degrees
        _, u, s = ow.dc3dwrapper(
            2.0 / 3.0,
            [0, xokada[i] + 5e3, yokada[i]],
            5e3,
            45,
            [-1e7, 1e7],
            [-10e3*sqrt(2) / 2, 10e3*sqrt(2) / 2],
            [0.0, 1.0, 0.0],
        )
        dispxokada[i] = u[2]
        dispyokada[i] = u[3]
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
        stressxxokada[i] = sxx
        stressyyokada[i] = syy
        stressxyokada[i] = sxy
    end

    # Off-fault displacements and stresses
    dispfaultconstvol, stressfaultconstvol = constdispstress(slip2dispstress, xokada, yokada, els, idx["fault"], faultslipconst[1:2:end], faultslipconst[2:2:end], mu, nu)
    dispfreesurfaceconstvol, stressfreesurfaceconstvol = constdispstress(slip2dispstress, xokada, yokada, els, idx["freesurf"], dispfreesurfaceconst[1:2:end], dispfreesurfaceconst[2:2:end], mu, nu)
    dispconst = dispfaultconstvol - dispfreesurfaceconstvol # Note negative sign
    stressconst = stressfaultconstvol - stressfreesurfaceconstvol # Note negative sign

    qux = transpose(reshape(dispfreesurfacequad[1:2:end], 3, nfreesurf))
    quy = transpose(reshape(dispfreesurfacequad[2:2:end], 3, nfreesurf))
    dispfaultquadvol, stressfaultquadvol = quaddispstress(slip2dispstress, xokada, yokada, els, idx["fault"], transpose(faultslipquad[1:2:end]), transpose(faultslipquad[2:2:end]), mu, nu)
    stressfreesurfacequadvol, stressfreesurfacequadvol = quaddispstress(slip2dispstress, xokada, yokada, els, idx["freesurf"], qux, quy, mu, nu)
    dispquad = dispfaultquadvol - dispfreesurfacequadvol # Note negative sign
    stressquad = stressfaultquadvol - stressfreesurfacequadvol # Note negative sign

    # Plot ux and uy profiles
    fontsize = 24
    markersize = 15
    linewidth = 2.0
    close("all")
    figure(figsize = (30, 15))

    ax = subplot(2, 3, 1)
    plot(xokada, log10.(abs.(stressconst[:, 1])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(stressquad[:, 1])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(stressxxokada[:, 1])), ":k", linewidth=2.0, label="Okada")
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([0, 8])
    gca().set_xticks([])
    gca().set_yticks([0, 4, 8])
    legend(fontsize=fontsize, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    # xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$\log \, \sigma_{xx}$ (Pa)", fontsize=fontsize)

    ax = subplot(2, 3, 2)
    plot(xokada, log10.(abs.(stressconst[:, 2])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(stressquad[:, 2])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(stressyyokada[:, 1])), ":k", linewidth=2.0, label="Okada")
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([0, 8])
    gca().set_xticks([])
    gca().set_yticks([0, 4, 8])
    legend(fontsize=fontsize, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$\log \, \sigma_{yy}$ (Pa)", fontsize=fontsize)

    ax = subplot(2, 3, 3)
    plot(xokada, log10.(abs.(stressconst[:, 3])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(stressquad[:, 3])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(stressxyokada[:, 1])), ":k", linewidth=2.0, label="Okada")
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([0, 8])
    gca().set_xticks([])
    gca().set_yticks([0, 4, 8])
    legend(fontsize=fontsize,frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$\log \, \sigma_{xy}$ (Pa)", fontsize=fontsize)

    ax = subplot(2, 3, 4)
    stressxxconsterror = 100 * (stressconst[:, 1] - stressxxokada[:, 1]) ./ stressxxokada[:, 1]
    stressxxquaderror = 100 * (stressquad[:, 1] - stressxxokada[:, 1]) ./ stressxxokada[:, 1]
    plot(xokada, log10.(abs.(stressxxconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM median %% error = %05.2f" median(abs.(σxxconsterror)))
    plot(xokada, log10.(abs.(stressxxquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM median %% error = %05.2f" median(abs.(σxxquaderror)))
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([-3, 6])
    gca().set_xticks([-50e3, 0, 50e3])
    gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsize, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$\log \, \sigma_{xx} \, \%$ error", fontsize=fontsize)

    ax = subplot(2, 3, 5)
    stressyyconsterror = 100 * (stressconst[:, 2] - stressyyokada[:, 1]) ./ stressyyokada[:, 1]
    stressyyquaderror = 100 * (stressquad[:, 2] - stressyyokada[:, 1]) ./ stressyyokada[:, 1]
    plot(xokada, log10.(abs.(stressyyconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM median %% error = %05.2f" median(abs.(σyyconsterror)))
    plot(xokada, log10.(abs.(stressyyquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM median %% error = %05.2f" median(abs.(σyyquaderror)))
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([-3, 6])
    gca().set_xticks([-50e3, 0, 50e3])
    gca().set_yticks([-3, 0, 3, 6])
    lh = legend(fontsize=fontsize, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{yy} \, \%$ error", fontsize=fontsize)

    ax = subplot(2, 3, 6)
    stressxyconsterror = 100 * (stressconst[:, 3] - stressxyokada[:, 1]) ./ stressxyokada[:, 1]
    stressxyquaderror = 100 * (stressquad[:, 3] - stressxyokada[:, 1]) ./ stressxyokada[:, 1]
    plot(xokada, log10.(abs.(stressxyconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM median %% error = %05.2f" median(abs.(σxyconsterror)))
    plot(xokada, log10.(abs.(stressxyquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM median %% error = %05.2f" median(abs.(σxyquaderror)))
    gca().set_xlim([-50e3, 50e3])
    gca().set_ylim([-3, 6])
    gca().set_xticks([-50e3, 0, 50e3])
    gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsize, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$\log \, \sigma_{xy} \, \%$ error", fontsize=fontsize)

    tight_layout()
    show()
end
ex_freesurface()
