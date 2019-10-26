using Revise
using PyCall
using PyPlot
using Bem2d


function plot6panel(els, xobs, yobs, npts, u, σ, titlelabel)
    σ = log10.(abs.(σ))

    fontsize = 20
    figure(figsize = (30, 20))
    ncontours = 50

    subplot(2, 3, 2)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(u[:, 1], npts, npts), ncontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$u_x$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(u[:, 1], npts, npts), ncontours, linewidths = 0.5, colors = "gray")
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$y$ (m)")

    subplot(2, 3, 3)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(u[:, 2], npts, npts), ncontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$u_y$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(u[:, 2], npts, npts), ncontours, linewidths = 0.5, colors = "gray")
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 3, 4)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(σ[:, 1], npts, npts), ncontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\sigma_{xx}$ (Pa)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(σ[:, 1], npts, npts), ncontours, linewidths = 0.5, colors = "gray")
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 3, 5)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(σ[:, 2], npts, npts), ncontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\sigma_{yy}$ (Pa)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(σ[:, 2], npts, npts), ncontours, linewidths = 0.5, colors = "gray")
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 3, 6)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(σ[:, 3], npts, npts), ncontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\sigma_{xy}$ (Pa)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(σ[:, 3], npts, npts), ncontours, linewidths = 0.5, colors = "gray")
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    suptitle(titlelabel, fontsize=2*fontsize)
    show()
end

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


    # Attempt at showing actual BEM slip
    faultidx = findall(x->x == "fault", els.name)
    freesurfidx = findall(x->x == "freesurf", els.name)
    uconstfault, σconstfault = constuσ(slip2uσ, xokada, zeros(size(xokada)), els, faultidx, faultslipconst[1:2:end], faultslipconst[2:2:end], μ, ν)
    uconstfreesurf, σconstfreesurf = constuσ(slip2uσ, xokada, zeros(size(xokada)), els, freesurfidx, ufreesurfaceconst[1:2:end], ufreesurfaceconst[2:2:end], μ, ν)
    uconst = uconstfault

    # Plot ux and uy profiles
    fontsize = 24
    markersize = 15
    linewidth = 0.5
    close("all")
    figure(figsize = (14, 18))

    ax = subplot(2, 1, 1)
    plot(xokada, uxokada, "-k", linewidth=linewidth, label="Okada")
    plot(xokada, uconst[1:2:end], "-g", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")

    plot(xplotconst, ufreesurfaceconst[1:2:end], "bo", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    plot(xplotquad, ufreesurfacequad[1:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "quad halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([-1.00, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    ylabel(L"$u_x$ (m)", fontsize=fontsize)

    ax = subplot(2, 1, 2)
    plot(xokada, uyokada, "-k", linewidth=linewidth, label="Okada")
    plot(xplotconst, ufreesurfaceconst[2:2:end], "bo", markeredgewidth=linewidth, markersize=markersize, label = "const halfspace")
    plot(xplotquad, ufreesurfacequad[2:2:end], "r+", markeredgewidth=linewidth, markersize=markersize, label = "quad halfspace")
    gca().set_xlim([-5, 5]); gca().set_ylim([-1.0, 1.0])
    gca().set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
    gca().set_yticks([-1.00, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00])
    legend(fontsize=fontsize)
    ax.tick_params("both", labelsize = fontsize)
    xlabel(L"$x$ (m)", fontsize=fontsize); ylabel(L"$u_y$ (m)", fontsize=fontsize)
    show()

    # # Now do volume solution to assess the effects of discontinuities
    # npts = 50
    # xobs, yobs = obsgrid(-5, -2, 5, 2, npts)
    #
    # faultidx = findall(x->x == "fault", els.name)
    # freesurfidx = findall(x->x == "freesurf", els.name)
    #
    # ufaultconstvol, σfaultconstvol = constuσ(slip2uσ, xobs, yobs, els, faultidx, faultslipconst[1:2:end], faultslipconst[2:2:end], μ, ν)
    # ufreesurfaceconstvol, σfreesurfaceconstvol = constuσ(slip2uσ, xobs, yobs, els, freesurfidx, ufreesurfaceconst[1:2:end], ufreesurfaceconst[2:2:end], μ, ν)
    # # plot6panel(els, xobs, yobs, npts, ufaultconstvol, σfaultconstvol, "fault only (CS elements)")
    # plot6panel(els, xobs, yobs, npts, ufreesurfaceconstvol, σfreesurfaceconstvol, "free surface only (CS elements)")
    # # plot6panel(els, xobs, yobs, npts, ufreesurfaceconstvol + ufaultconstvol, σfreesurfaceconstvol + σfaultconstvol, "total (CS elements)")
    #
    # qux = transpose(reshape(ufreesurfacequad[1:2:end], 3, nfreesurf))
    # quy = transpose(reshape(ufreesurfacequad[2:2:end], 3, nfreesurf))
    #
    # ufaultquadvol, σfaultquadvol = quaduσ(slip2uσ, xobs, yobs, els, faultidx, transpose(faultslipquad[1:2:end]), transpose(faultslipquad[2:2:end]), μ, ν)
    # ufreesurfacequadvol, σfreesurfacequadvol = quaduσ(slip2uσ, xobs, yobs, els, freesurfidx, qux, quy, μ, ν)
    # # plot6panel(els, xobs, yobs, npts, ufaultquadvol, σfaultquadvol, "fault only (3QN elements)")
    # plot6panel(els, xobs, yobs, npts, ufreesurfacequadvol, σfreesurfacequadvol, "free surface only (3QN elements)")
    # # plot6panel(els, xobs, yobs, npts, ufreesurfacequadvol + ufaultquadvol, σfreesurfacequadvol + σfaultquadvol, "total (3QN elements)")


end
ex_freesurface()
