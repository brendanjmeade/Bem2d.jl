using Revise
using PyCall
using PyPlot
using Bem2d

function ex_freesurface()
    mu = 30e9
    nu = 0.25

    # Free surface
    els = Elements(Int(1e5))
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

    # Conviencience structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    partialsquad = initpartials(els)

    # Constant slip fault
    partialsconst["disp"]["fault"]["freesurf"], partialsconst["stress"]["fault"]["freesurf"], partialsconst["trac"]["fault"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    partialsconst["disp"]["freesurf"]["freesurf"], partialsconst["stress"]["freesurf"]["freesurf"], partialsconst["trac"]["freesurf"]["freesurf"] = partialsconstdispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    faultslipconst = sqrt(2) / 2 * [1 ; 1]
    ufullspaceconst = partialsconst["disp"]["fault"]["freesurf"] * faultslipconst
    ufreesurfaceconst = inv(partialsconst["trac"]["freesurf"]["freesurf"]) * (partialsconst["trac"]["fault"]["freesurf"] * faultslipconst)
    xplotconst = els.xcenter[idx["freesurf"]]

    # Quadratic slip fault
    partialsquad["disp"]["fault"]["freesurf"], partialsquad["stress"]["fault"]["freesurf"], partialsquad["trac"]["fault"]["freesurf"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["freesurf"], mu, nu)
    partialsquad["disp"]["freesurf"]["freesurf"], partialsquad["stress"]["freesurf"]["freesurf"], partialsquad["trac"]["freesurf"]["freesurf"] = partialsquaddispstress(slip2dispstress, els, idx["freesurf"], idx["freesurf"], mu, nu)
    xplotquad = sort(els.xnodes[idx["freesurf"], :][:])
    faultslipquad = sqrt(2) ./ 2 * ones(6)
    ufullspacequad = partialsquad["disp"]["fault"]["freesurf"] * faultslipquad
    ufreesurfacequad = inv(partialsquad["trac"]["freesurf"]["freesurf"]) * (partialsquad["trac"]["fault"]["freesurf"] * faultslipquad)

    # Okada solution
    ow = pyimport("okada_wrapper") # from okada_wrapper import dc3dwrapper
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