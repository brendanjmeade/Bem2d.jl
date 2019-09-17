using Revise
using PyCall
using PyPlot
using ColorSchemes
using Bem2d

function plotfunction(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad, xobs, yobs, uconstfarfield, σconstfarfield)
    cmap = ColorScheme([Colors.RGB(255.0 / 255.0, 20.0 / 255.0, 40.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 47.0 / 255.0, 146.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 138.0 / 255.0, 216.0 / 255.0),
                        Colors.RGB(55.0 / 255.0, 145.0 / 255.0, 230.0 / 255.0),
                        Colors.RGB(150.0 / 255.0, 230.0 / 255.0, 80.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 251.0 / 255.0, 0.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 255.0 / 255.0, 255.0 / 255.0)],
                        "rycroft", "psychedelic")
    cmap = ColorScheme([get(cmap, i) for i in 0.0:0.001:1.0])

    scalestress = 1e6
    vmin = -1e6 / scalestress
    vmax = 1e6 / scalestress
    ncontours = 50
    fontsize = 6
    xmin = -5e3
    xmax = 5e3

    figure(figsize = (10, 6))
    ax = subplot(2, 3, 1)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, vmin = vmin, vmax = vmax, cmap = ColorMap(cmap.colors))
    clim(vmin, vmax)
    cbar = colorbar(orientation = "horizontal", fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize) 
    contour(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, vmin = vmin, vmax = vmax,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize = fontsize); xlim([xmin, xmax]); xticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(2, 3, 2)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 2] ./ scalestress, size(xobs)), ncontours, vmin = vmin, vmax = vmax, cmap = ColorMap(cmap.colors))
    clim(vmin, vmax)
    cbar = colorbar(orientation = "horizontal", fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize) 
    contour(xobs, yobs, reshape(σconstfarfield[:, 2] ./ scalestress, size(xobs)), ncontours, vmin = vmin, vmax = vmax,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(2, 3, 3)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 3] ./ scalestress, size(xobs)), ncontours, vmin = vmin, vmax = vmax, cmap = ColorMap(cmap.colors))
    clim(vmin, vmax)
    cbar = colorbar(orientation = "horizontal", fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize) 
    contour(xobs, yobs, reshape(σconstfarfield[:, 3] ./ scalestress, size(xobs)), ncontours, vmin = vmin, vmax = vmax,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(2, 3, 4)
    plot(xcenters, σconst[1:3:end] ./ scalestress, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{xx}$ constant", zorder = 2)
    plot(xnodes, σquad[1:3:end] ./ scalestress, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{xx}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress])

    ax = subplot(2, 3, 5)
    plot(xcenters, σconst[2:3:end] ./ scalestress, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{yy}$ constant", zorder = 2)
    plot(xnodes, σquad[2:3:end] ./ scalestress, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{yy}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize);
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress]); yticks([])

    ax = subplot(2, 3, 6)
    plot(xcenters, σconst[3:3:end] ./ scalestress, ".b", markerfacecolor = "navy", markeredgewidth = 0.0, label = L"$σ_{xy}$ constant", zorder = 2)
    plot(xnodes, σquad[3:3:end] ./ scalestress, "or", markerfacecolor = "coral", markeredgewidth = 0.0, label = L"$σ_{xy}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize);    # ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress]); yticks([])
    tight_layout(); show()
end

function fig_constquadstresscomp()
    # Material and geometric constants
    μ = 3e10
    ν = 0.25
    nels = 20
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
    npts = 100; obswidth = 5e3
    xobs, yobs = obsgrid(-obswidth, 0, obswidth, 5e3, npts)

    # Evaluation points
    xcenters = els.xcenter[1:els.endidx]
    ycenters = els.ycenter[1:els.endidx]
    xnodes = sort(els.xnodes[1:els.endidx, :][:])
    ynodes = sort(els.ynodes[1:els.endidx, :][:])

    # Slip for volumetric models
    slope = 0.0001
    constxslip = zeros(2, nels)
    constyslip = zeros(2, nels)
    constxslip[1, :] .= 1.0
    constxslip[2, :] = slope .* xcenters

    # # Linear x-slip only
    # slope = 0.001
    # constxslip = slope .* xcenters
    # constyslip = zeros(nels)
    # quadxslip = transpose(reshape(slope .* xnodes, 3, nels))
    # quadyslip = zeros(size(quadxslip))


    # quadxslip = zeros(2, nels, 3)
    # quadyslip = zeros(2, nels, 3)
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
    for i in 1:size(slipconst)[1]
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
