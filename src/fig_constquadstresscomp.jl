using Revise
using PyCall
using PyPlot
using Bem2d

function plotfunction(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad, xobs, yobs, uconstfarfield, σconstfarfield)
    scalestress = 1e6
    ncontours = 30
    fontsize = 6
    xmin = -5e3
    xmax = 5e3
    cmap = rycroftcmap()
    
    figure(figsize = (10, 6))
    ax = subplot(2, 3, 1)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize) 
    contour(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize = fontsize); xlim([xmin, xmax]); xticks([])
    yticks([-2500, 0, 2500]);
    ax.tick_params(labelsize = fontsize)

    ax = subplot(2, 3, 2)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize) 
    contour(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(2, 3, 3)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 3] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize) 
    contour(xobs, yobs, reshape(σconstfarfield[:, 3] ./ scalestress, size(xobs)), ncontours,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(2, 3, 4)
    plot(xcenters, σconst[1:3:end] ./ scalestress, ".k", markerfacecolor = "k", markeredgewidth = 0.0, label = L"$σ_{xx}$ constant", zorder = 2)
    plot(xnodes, σquad[1:3:end] ./ scalestress, "ok", markerfacecolor = "lightgray", markeredgewidth = 0.25, label = L"$σ_{xx}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize); ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    yticks([-10, 0, 10]);
    ylim([-1e7 / scalestress, 1e7 / scalestress])

    ax = subplot(2, 3, 5)
    plot(xcenters, σconst[2:3:end] ./ scalestress, ".k", markerfacecolor = "k", markeredgewidth = 0.0, label = L"$σ_{yy}$ constant", zorder = 2)
    plot(xnodes, σquad[2:3:end] ./ scalestress, "ok", markerfacecolor = "lightgray", markeredgewidth = 0.25, label = L"$σ_{yy}$ quadratic", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize); xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize);
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress]); yticks([])

    ax = subplot(2, 3, 6)
    plot(xcenters, σconst[3:3:end] ./ scalestress, ".k", markerfacecolor = "k", markeredgewidth = 0.0, label = L"$σ_{xy}$ constant", zorder = 2)
    plot(xnodes, σquad[3:3:end] ./ scalestress, "ok", markerfacecolor = "lightgray", markeredgewidth = 0.25, label = L"$σ_{xy}$ quadratic", zorder = 1)
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
    nels = 10
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
    xobs, yobs = obsgrid(-obswidth, -2.5e3, obswidth, 2.5e3, npts)

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
