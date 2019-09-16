using Revise
using PyCall
using PyPlot
using ColorSchemes
using Bem2d

function plotfunction(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad, xobs, yobs, uconstfarfield, σconstfarfield)
    scalestress = 1e6
    cmap = ColorScheme([Colors.RGB(255.0 / 255.0, 20.0 / 255.0, 40.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 47.0 / 255.0, 146.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 138.0 / 255.0, 216.0 / 255.0),
                        Colors.RGB(55.0 / 255.0, 145.0 / 255.0, 230.0 / 255.0),
                        Colors.RGB(150.0 / 255.0, 230.0 / 255.0, 80.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 251.0 / 255.0, 0.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 255.0 / 255.0, 255.0 / 255.0)],
                        "rycroft", "psychedelic")
    cmap = ColorScheme([get(cmap, i) for i in 0.0:0.001:1.0])
    # contourf(rand(30, 30), 50, cmap = ColorMap(cmap.colors))
    # show()
    fontsize = 6
    figure(figsize = (6, 6))
    ax = subplot(3, 3, 1)

    xlimit = [minimum(xobs) maximum(xobs)]
    ylimit = [minimum(yobs) maximum(yobs)]
    scale = 5e-1
    ncontours = 50
    # field = σconstfarfield[:, 1]
    # fieldmax = maximum(@.abs(field))
    contourf(xobs, yobs, reshape(σconstfarfield[:, 1], size(xobs)), ncontours, vmin = -scale * fieldmax, vmax = scale * fieldmax, cmap = ColorMap(cmap.colors))
    # clim(-scale * fieldmax, scale * fieldmax)
    cbar = colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize) 
    contour(xobs, yobs, reshape(σconstfarfield[:, 1], size(xobs)), ncontours, vmin = -scale * fieldmax, vmax = scale * fieldmax, linewidths = 0.25, colors = "w", linestyles = "-")
    gca().set_aspect("equal")
    # for i in 1:els.endidx
    #     plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-b", color = "b", linewidth = 0.5)
    # end
    ylabel("y (m)", fontsize = fontsize); xlim([-20000, 20000]); xticks([-20000, 0, 20000])
    ax.tick_params(labelsize = fontsize)

    # ax = subplot(3, 3, 2)
    # plot(xcenters, uconst[1:2:end], ".b",  markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$u_x$ constant", zorder = 2)
    # plot(xnodes, uquad[1:2:end], "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$u_x$ quadratic", zorder = 1)
    # legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    # ylabel(L"$u$ (m)", fontsize = fontsize)
    # ax[:tick_params]("both", labelsize = fontsize)
    # ylim([-0.6, 0.6])

    # ax = subplot(3, 3, 3)
    # plot(xcenters, uconst[2:2:end], ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$u_y$ constant", zorder = 2)
    # plot(xnodes, uquad[2:2:end], "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$u_y$ quadratic", zorder = 1)
    # legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    # ylabel(L"$u$ (m)", fontsize = fontsize)
    # ax[:tick_params]("both", labelsize = fontsize)
    # ylim([-0.6, 0.6])

    ax = subplot(3, 3, 7)
    plot(xcenters, σconst[1:3:end] ./ scalestress, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{xx}$ constant", zorder = 2)
    plot(xnodes, σquad[1:3:end] ./ scalestress, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{xx}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress])

    ax = subplot(3, 3, 8)
    plot(xcenters, σconst[2:3:end] ./ scalestress, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{yy}$ constant", zorder = 2)
    plot(xnodes, σquad[2:3:end] ./ scalestress, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{yy}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress])

    ax = subplot(3, 3, 9)
    plot(xcenters, σconst[3:3:end] ./ scalestress, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{xy}$ constant", zorder = 2)
    plot(xnodes, σquad[3:3:end] ./ scalestress, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{xy}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([-10000, 10000]); xticks([-10000, 0, 10000])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress])
    tight_layout(); show()
end

function fig_constquadstresscomp()
    # Material and geometric constants
    μ = 3e10
    ν = 0.25
    nels = 5
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

    # Observation coordinates for plotting
    npts = 100; obswidth = 20e3
    xobs, yobs = obsgrid(-obswidth, 0, obswidth, obswidth, npts)
    constxslip = zeros(2, nels)
    constyslip = zeros(2, nels)
    constxslip[1, :] .= 1.0
    # quadxslip = zeros(2, 3 * nels)
    # quadyslip = zeros(2, 3 * nels)
    # quadxslip = repeat(constxslip, 1, 3)
    # quadyslip = repeat(constyslip, 1, 3)

    # Evaluation points and slip
    xcenters = els.xcenter[1:els.endidx]
    ycenters = els.ycenter[1:els.endidx]
    xnodes = sort(els.xnodes[1:els.endidx, :][:])
    ynodes = sort(els.ynodes[1:els.endidx, :][:])

    # Parameters for the test cases
    slipconst = zeros(2, 2 * nels)
    slipquad = zeros(2, 6 * nels)

    # Constant strike-slip
    slipconst[1, 1:2:end] .= 1
    slipquad[1, 1:2:end] .= 1

    # Linear strike-slip
    slope = 0.0001
    slipconst[2, 1:2:end] = slope * xcenters
    slipquad[2, 1:2:end] = slope * xnodes

    # Predict on-fault displacements, stresses, and tractions
    close("all")
    for i in 1:1#size(slipconst)[1]
        # On fault evaluation
        uconst = ∂uconst * slipconst[i, :]
        σconst = ∂σconst * slipconst[i, :]
        tconst = ∂tconst * slipconst[i, :]
        uquad = ∂uquad * slipquad[i, :]
        σquad = ∂σquad * slipquad[i, :]
        tquad = ∂tquad * slipquad[i, :]

        # Far-field evaluation
        uconstfarfield, σconstfarfield = constuσ(slip2uσ, xobs, yobs, els, srcidx, constxslip[i, :], constyslip[i, :], μ, ν)
        # uquadfarfield, σquadfarfield = quaduσ(slip2uσ, xobs, yobs, els, srcidx, quadxslip[i, :], quadyslip[i, :], μ, ν)
    
        plotfunction(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad, reshape(xobs, npts, npts), reshape(yobs, npts, npts), uconstfarfield, σconstfarfield)
    end

end
fig_constquadstresscomp()
