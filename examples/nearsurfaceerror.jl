using Revise
using Printf
using PyCall
using PyPlot
using Statistics
using Bem2d

"""
    nearsurfaceerror()

Compare BEM errors from constant and quadratic elements.  The reference
model is an Okada approximation of a 2d thrust fault.  BEM errors are
evaluated not on any boundary but rather at a small distance away from
the free surface boundary.
"""
function nearsurfaceerror()
    mu = 30e9
    nu = 0.25

    # Free surface
    els = Elements(Int(1e5))
    # nfreesurf = 20
    # x1, y1, x2, y2 = discretizedline(-50e3, 0, 50e3, 0, nfreesurf)
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-100e3, 0, 100e3, 0, nfreesurf)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "freesurf"
    end
    standardize_elements!(els)

    @show "element width"
    @show x2[1] - x1[1]


    x1, y1, x2, y2 = discretizedline(-100e3, 0, 100e3, 0, 3 * nfreesurf)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "freesurf2"
    end
    standardize_elements!(els)

    
    # 45 degree dipping fault
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault)
    y1 = y1 .- 5e3
    y2 = y2 .- 5e3
    
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
    partialsconst2 = initpartials(els)
    partialsquad = initpartials(els)

    # Constant slip fault
    partialsconst["disp"]["fault"]["freesurf"], _, partialsconst["trac"]["fault"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    partialsconst["disp"]["freesurf"]["freesurf"], _, partialsconst["trac"]["freesurf"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    faultslipconst = sqrt(2) / 2 * [1 ; 1]
    dispfullspaceconst = partialsconst["disp"]["fault"]["freesurf"] * faultslipconst
    dispfreesurfaceconst = inv(partialsconst["trac"]["freesurf"]["freesurf"]) * (partialsconst["trac"]["fault"]["freesurf"] * faultslipconst)
    xplotconst = els.xcenter[idx["freesurf"]]

    # HIGH RES: Constant slip fault
    partialsconst2["disp"]["fault"]["freesurf2"], _, partialsconst2["trac"]["fault"]["freesurf2"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurf2"], mu, nu)
    partialsconst2["disp"]["freesurf2"]["freesurf2"], _, partialsconst2["trac"]["freesurf2"]["freesurf2"] = partialsconstdispstress(slip2dispstress, els, idx["freesurf2"], idx["freesurf2"], mu, nu)
    faultslipconst2 = sqrt(2) / 2 * [1 ; 1]
    dispfullspaceconst2 = partialsconst2["disp"]["fault"]["freesurf2"] * faultslipconst2
    dispfreesurfaceconst2 = inv(partialsconst2["trac"]["freesurf2"]["freesurf2"]) * (partialsconst2["trac"]["fault"]["freesurf2"] * faultslipconst2)
    xplotconst2 = els.xcenter[idx["freesurf2"]]

    # Quadratic slip fault
    partialsquad["disp"]["fault"]["freesurf"], _, partialsquad["trac"]["fault"]["freesurf"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    partialsquad["disp"]["freesurf"]["freesurf"], _, partialsquad["trac"]["freesurf"]["freesurf"] = partialsquaddispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    xplotquad = sort(els.xnodes[idx["freesurf"], :][:])
    faultslipquad = sqrt(2) ./ 2 * ones(6)
    dispfullspacequad = partialsquad["disp"]["fault"]["freesurf"] * faultslipquad
    dispfreesurfacequad = inv(partialsquad["trac"]["freesurf"]["freesurf"]) * (partialsquad["trac"]["fault"]["freesurf"] * faultslipquad)
    
    # Okada solution
    ow = pyimport("okada_wrapper")# from okada_wrapper import dc3dwrapper
    xokada = collect(LinRange(-2e3, 2e3, 10000))
    offset = abs(els.xnodes[1, 2] - els.xnodes[1, 1])
    offset = 100
    @show offset
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
            # 5e3,
            10e3,            
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
        sxx = mu * (exx + eyy) + 2 * mu * exx
        syy = mu * (exx + eyy) + 2 * mu * eyy
        sxy = 2 * mu * exy
        stressxxokada[i] = sxx
        stressyyokada[i] = syy
        stressxyokada[i] = sxy
    end

    # Off-fault displacements and stresses
    dispfaultconstvol, stressfaultconstvol = constdispstress(slip2dispstress, xokada, yokada, els, idx["fault"], faultslipconst[1:2:end], faultslipconst[2:2:end], mu, nu)
    dispfreesurfaceconstvol, stressfreesurfaceconstvol = constdispstress(slip2dispstress, xokada, yokada, els, idx["freesurf"], dispfreesurfaceconst[1:2:end], dispfreesurfaceconst[2:2:end], mu, nu)
    dispconst = dispfaultconstvol - dispfreesurfaceconstvol # Note negative sign
    stressconst = stressfaultconstvol - stressfreesurfaceconstvol # Note negative sign

    # HIGH RES: Off-fault displacements and stresses
    dispfaultconstvol2, stressfaultconstvol2 = constdispstress(slip2dispstress, xokada, yokada, els, idx["fault"], faultslipconst[1:2:end], faultslipconst[2:2:end], mu, nu)
    dispfreesurfaceconstvol2, stressfreesurfaceconstvol2 = constdispstress(slip2dispstress, xokada, yokada, els, idx["freesurf2"], dispfreesurfaceconst2[1:2:end], dispfreesurfaceconst2[2:2:end], mu, nu)
    dispconst2 = dispfaultconstvol2 - dispfreesurfaceconstvol2 # Note negative sign
    stressconst2 = stressfaultconstvol2 - stressfreesurfaceconstvol2 # Note negative sign

    qux = transpose(reshape(dispfreesurfacequad[1:2:end], 3, nfreesurf))
    quy = transpose(reshape(dispfreesurfacequad[2:2:end], 3, nfreesurf))
    dispfaultquadvol, stressfaultquadvol = quaddispstress(slip2dispstress, xokada, yokada, els, idx["fault"], transpose(faultslipquad[1:2:end]), transpose(faultslipquad[2:2:end]), mu, nu)
    dispfreesurfacequadvol, stressfreesurfacequadvol = quaddispstress(slip2dispstress, xokada, yokada, els, idx["freesurf"], qux, quy, mu, nu)
    dispquad = dispfaultquadvol - dispfreesurfacequadvol # Note negative sign
    stressquad = stressfaultquadvol - stressfreesurfacequadvol # Note negative sign

    # Plot ux and uy profiles
    xmin = -2e3
    xmax = 2e3
    fontsize = 24
    fontsizelegend = 20
    markersize = 15
    linewidth = 2.0
    close("all")
    figure(figsize = (30, 20))

    ax = subplot(3, 3, 1)
    plot(xokada, log10.(abs.(stressconst[:, 1])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(stressconst2[:, 1])), "-m", linewidth=linewidth, label="3 x CS BEM")
    plot(xokada, log10.(abs.(stressquad[:, 1])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(stressxxokada)), "--k", linewidth=3.0, label="analytic")
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([2, 8])
    gca().set_xticks([xmin, 0, xmax])
    # gca().set_yticks([0, 4, 8])
    legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$\log \, |\sigma_{xx}|$ (Pa)", fontsize=fontsize)

    ax = subplot(3, 3, 2)
    plot(xokada, log10.(abs.(stressconst[:, 2])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(stressconst2[:, 2])), "-m", linewidth=linewidth, label="3 x CS BEM")
    plot(xokada, log10.(abs.(stressquad[:, 2])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(stressyyokada[:, 1])), "--k", linewidth=3.0, label="analytic")
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([2, 8])
    gca().set_xticks([xmin, 0, xmax])
    # gca().set_yticks([0, 4, 8])
    legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$\log \, |\sigma_{yy}|$ (Pa)", fontsize=fontsize)

    ax = subplot(3, 3, 3)
    plot(xokada, log10.(abs.(stressconst[:, 3])), "-c", linewidth=linewidth, label="CS BEM")
    plot(xokada, log10.(abs.(stressconst2[:, 3])), "-m", linewidth=linewidth, label="3 x CS BEM")
    plot(xokada, log10.(abs.(stressquad[:, 3])), "-r", linewidth=linewidth, label="3QN BEM")
    plot(xokada, log10.(abs.(stressxyokada[:, 1])), "--k", linewidth=3.0, label="analytic")
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([2, 8])
    gca().set_xticks([xmin, 0, xmax])
    # gca().set_yticks([0, 4, 8])
    legend(fontsize=fontsizelegend,frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$\log \, |\sigma_{xy}|$ (Pa)", fontsize=fontsize)

    # Absolute error
    ax = subplot(3, 3, 4)
    stressxxconsterror = abs.(stressconst[:, 1] - stressxxokada[:, 1])
    stressxxconsterror2 = abs.(stressconst2[:, 1] - stressxxokada[:, 1])
    stressxxquaderror = abs.(stressquad[:, 1] - stressxxokada[:, 1])
    plot(xokada, log10.(abs.(stressxxconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM mean error = %05.2f" mean(abs.(stressxxconsterror)))
    plot(xokada, log10.(abs.(stressxxconsterror2)), "-m", linewidth=linewidth, label=@sprintf "3 x CS BEM mean error = %05.2f" mean(abs.(stressxxconsterror2)))
    plot(xokada, log10.(abs.(stressxxquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM mean error = %05.2f" mean(abs.(stressxxquaderror)))
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([2, 8])
    gca().set_xticks([xmin, 0, xmax])
    # gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$|\sigma_{xx}$ error|", fontsize=fontsize)

    ax = subplot(3, 3, 5)
    stressyyconsterror = abs.(stressconst[:, 2] - stressyyokada[:, 1])
    stressyyconsterror2 = abs.(stressconst2[:, 2] - stressyyokada[:, 1])
    stressyyquaderror = abs.(stressquad[:, 2] - stressyyokada[:, 1])
    plot(xokada, log10.(abs.(stressyyconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM mean error = %05.2f" mean(abs.(stressyyconsterror)))
    plot(xokada, log10.(abs.(stressyyconsterror2)), "-m", linewidth=linewidth, label=@sprintf "3 x CS BEM mean error = %05.2f" mean(abs.(stressyyconsterror2)))
    plot(xokada, log10.(abs.(stressyyquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM mean error = %05.2f" mean(abs.(stressyyquaderror)))
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([2, 8])
    gca().set_xticks([xmin, 0, xmax])
    # gca().set_yticks([-3, 0, 3, 6])
    lh = legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$|\sigma_{yy}$ error|", fontsize=fontsize)

    ax = subplot(3, 3, 6)
    stressxyconsterror = abs.(stressconst[:, 3] - stressxyokada[:, 1])
    stressxyconsterror2 = abs.(stressconst2[:, 3] - stressxyokada[:, 1])
    stressxyquaderror = abs.(stressquad[:, 3] - stressxyokada[:, 1])
    plot(xokada, log10.(abs.(stressxyconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM mean error = %05.2f" mean(abs.(stressxyconsterror)))
    plot(xokada, log10.(abs.(stressxyconsterror2)), "-m", linewidth=linewidth, label=@sprintf "3 x CS BEM mean error = %05.2f" mean(abs.(stressxyconsterror2)))
    plot(xokada, log10.(abs.(stressxyquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM mean error = %05.2f" mean(abs.(stressxyquaderror)))
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([2, 8])
    gca().set_xticks([xmin, 0, xmax])
    # gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$|\sigma_{xy}$ error|", fontsize=fontsize)

    ax = subplot(3, 3, 7)
    stressxxconsterror = 100 * (stressconst[:, 1] - stressxxokada[:, 1]) ./ stressxxokada[:, 1]
    stressxxconsterror2 = 100 * (stressconst2[:, 1] - stressxxokada[:, 1]) ./ stressxxokada[:, 1]
    stressxxquaderror = 100 * (stressquad[:, 1] - stressxxokada[:, 1]) ./ stressxxokada[:, 1]
    plot(xokada, log10.(abs.(stressxxconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM median %% error = %05.2f" median(abs.(stressxxconsterror)))
    plot(xokada, log10.(abs.(stressxxconsterror2)), "-m", linewidth=linewidth, label=@sprintf "3 x CS BEM median %% error = %05.2f" median(abs.(stressxxconsterror2)))
    plot(xokada, log10.(abs.(stressxxquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM median %% error = %05.2f" median(abs.(stressxxquaderror)))
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([-3, 6])
    gca().set_xticks([xmin, 0, xmax])
    gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$\log \, |\sigma_{xx} \, \%$ error|", fontsize=fontsize)

    ax = subplot(3, 3, 8)
    stressyyconsterror = 100 * (stressconst[:, 2] - stressyyokada[:, 1]) ./ stressyyokada[:, 1]
    stressyyconsterror2 = 100 * (stressconst2[:, 2] - stressyyokada[:, 1]) ./ stressyyokada[:, 1]
    stressyyquaderror = 100 * (stressquad[:, 2] - stressyyokada[:, 1]) ./ stressyyokada[:, 1]
    plot(xokada, log10.(abs.(stressyyconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM median %% error = %05.2f" median(abs.(stressyyconsterror)))
    plot(xokada, log10.(abs.(stressyyconsterror2)), "-m", linewidth=linewidth, label=@sprintf "3 x CS BEM median %% error = %05.2f" median(abs.(stressyyconsterror2)))
    plot(xokada, log10.(abs.(stressyyquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM median %% error = %05.2f" median(abs.(stressyyquaderror)))
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([-3, 6])
    gca().set_xticks([xmin, 0, xmax])
    gca().set_yticks([-3, 0, 3, 6])
    lh = legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize);
    ylabel(L"$\log \, |\sigma_{yy} \, \%$ error|", fontsize=fontsize)

    ax = subplot(3, 3, 9)
    stressxyconsterror = 100 * (stressconst[:, 3] - stressxyokada[:, 1]) ./ stressxyokada[:, 1]
    stressxyconsterror2 = 100 * (stressconst2[:, 3] - stressxyokada[:, 1]) ./ stressxyokada[:, 1]
    stressxyquaderror = 100 * (stressquad[:, 3] - stressxyokada[:, 1]) ./ stressxyokada[:, 1]
    plot(xokada, log10.(abs.(stressxyconsterror)), "-c", linewidth=linewidth, label=@sprintf "CS BEM median %% error = %05.2f" median(abs.(stressxyconsterror)))
    plot(xokada, log10.(abs.(stressxyconsterror2)), "-m", linewidth=linewidth, label=@sprintf "3 x CS BEM median %% error = %05.2f" median(abs.(stressxyconsterror2)))
    plot(xokada, log10.(abs.(stressxyquaderror)), "-r", linewidth=linewidth, label=@sprintf "3QN BEM median %% error = %05.2f" median(abs.(stressxyquaderror)))
    gca().set_xlim([xmin, xmax])
    gca().set_ylim([-3, 6])
    gca().set_xticks([xmin, 0, xmax])
    gca().set_yticks([-3, 0, 3, 6])
    legend(fontsize=fontsizelegend, frameon=true, facecolor="white", framealpha=1.0, loc=1)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize);
    ylabel(L"$\log \, |\sigma_{xy} \, \%$ error|", fontsize=fontsize)

    tight_layout(pad=5.0)

    show()
end
nearsurfaceerror()
