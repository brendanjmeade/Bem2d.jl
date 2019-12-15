using Revise
using PyCall
using PyPlot
using Bem2d


function plotfunction(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad, xobs, yobs, uconstfarfield, σconstfarfield, uquadfarfield, σquadfarfield)
    scalestress = 1e6
    contourlevels = 5
    contourlevels = collect(LinRange(4, 7, 6))
    fontsize = 7
    xmin = -5e3
    xmax = 5e3
    cmap = rycroftcmap()

    figure(figsize = (6, 5), dpi=300)
    ax = subplot(3, 3, 1)
    mat = reshape(log10.(abs.(σconstfarfield[:, 1])), size(xobs))
    mat[findall(x -> x < minimum(contourlevels), mat)] .= minimum(contourlevels)
    mat[findall(x -> x > maximum(contourlevels), mat)] .= maximum(contourlevels)
    contourf(xobs, yobs, mat, contourlevels, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\log\|\sigma_{xx}\|$ (Pa) CS", fontsize = fontsize)
    contour(xobs, yobs, mat, contourlevels, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize = fontsize); xlim([xmin, xmax]); xticks([])
    yticks([-2500, 0, 2500]);
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 2)
    mat = reshape(log10.(abs.(σconstfarfield[:, 2])), size(xobs))
    mat[findall(x -> x < minimum(contourlevels), mat)] .= minimum(contourlevels)
    mat[findall(x -> x > maximum(contourlevels), mat)] .= maximum(contourlevels)
    contourf(xobs, yobs, mat, contourlevels, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\log\|\sigma_{yy}\|$ (Pa) CS", fontsize = fontsize)
    contour(xobs, yobs, mat, contourlevels, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 3)
    mat = reshape(log10.(abs.(σconstfarfield[:, 3])), size(xobs))
    mat[findall(x -> x < minimum(contourlevels), mat)] .= minimum(contourlevels)
    mat[findall(x -> x > maximum(contourlevels), mat)] .= maximum(contourlevels)
    contourf(xobs, yobs, mat, contourlevels, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\log\|\sigma_{xy}\|$ (Pa) CS", fontsize = fontsize)
    contour(xobs, yobs, mat, contourlevels, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 4)
    mat = reshape(log10.(abs.(σquadfarfield[:, 1])), size(xobs))
    mat[findall(x -> x < minimum(contourlevels), mat)] .= minimum(contourlevels)
    mat[findall(x -> x > maximum(contourlevels), mat)] .= maximum(contourlevels)
    contourf(xobs, yobs, mat, contourlevels, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\log\|\sigma_{xx}\|$ (Pa) 3QN", fontsize = fontsize)
    contour(xobs, yobs, mat, contourlevels, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize = fontsize); xlim([xmin, xmax]); xticks([])
    yticks([-2500, 0, 2500]);
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 5)
    mat = reshape(log10.(abs.(σquadfarfield[:, 2])), size(xobs))
    mat[findall(x -> x < minimum(contourlevels), mat)] .= minimum(contourlevels)
    mat[findall(x -> x > maximum(contourlevels), mat)] .= maximum(contourlevels)
    contourf(xobs, yobs, mat, contourlevels, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\log\|\sigma_{yy}\|$ (Pa) 3QN", fontsize = fontsize)
    contour(xobs, yobs, mat, contourlevels, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 6)
    mat = reshape(log10.(abs.(σquadfarfield[:, 3])), size(xobs))
    mat[findall(x -> x < minimum(contourlevels), mat)] .= minimum(contourlevels)
    mat[findall(x -> x > maximum(contourlevels), mat)] .= maximum(contourlevels)
    contourf(xobs, yobs, mat, contourlevels, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\log\|\sigma_{xy}\|$ (Pa) 3QN", fontsize = fontsize)
    contour(xobs, yobs, mat, contourlevels,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 7)
    plot(xcenters, σconst[1:3:end] ./ scalestress, ".k", markerfacecolor = "blue", markeredgewidth = 0.0, label = L"$σ_{xx}$ CS", zorder = 2)
    plot(xnodes, σquad[1:3:end] ./ scalestress, "ok", markerfacecolor = "red", markeredgewidth = 0.25, label = L"$σ_{xx}$ 3QN", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize, frameon=false)
    xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize)
    ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    yticks([-10, 0, 10])
    ylim([-1e7 / scalestress, 1e7 / scalestress])

    ax = subplot(3, 3, 8)
    plot(xcenters, σconst[2:3:end] ./ scalestress, ".k", markerfacecolor = "blue", markeredgewidth = 0.0, label = L"$σ_{yy}$ CS", zorder = 2)
    plot(xnodes, σquad[2:3:end] ./ scalestress, "ok", markerfacecolor = "red", markeredgewidth = 0.25, label = L"$σ_{yy}$ 3QN", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize, frameon=false)
    xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress]); yticks([])

    ax = subplot(3, 3, 9)
    plot(xcenters, σconst[3:3:end] ./ scalestress, ".k", markerfacecolor = "blue", markeredgewidth = 0.0, label = L"$σ_{xy}$ CS", zorder = 2)
    plot(xnodes, σquad[3:3:end] ./ scalestress, "ok", markerfacecolor = "red", markeredgewidth = 0.25, label = L"$σ_{xy}$ 3QN", zorder = 1)
    legend(loc = "lower right", fontsize = fontsize, frameon=false)
    xlim([xmin, xmax]); xticks([xmin, 0, xmax])
    xlabel(L"$x$ (m)", fontsize = fontsize)    # ylabel(L"$\sigma$ (MPa)", fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ylim([-1e7 / scalestress, 1e7 / scalestress]); yticks([])
    tight_layout()
    show()
end



function fig_constquadstresscomp()
    # Material and geometric constants
    mu = 3e10
    nu = 0.25
    nels = 20
    els = Elements(Int(1e5))
    L = 10000
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nels)
    x2[11] = 500.0
    x1[12] = 500.0
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Convenience data structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    partialsquad = initpartials(els)

    partialsconst["disp"]["fault"]["fault"], partialsconst["stress"]["fault"]["fault"], partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    partialsquad["disp"]["fault"]["fault"], partialsquad["stress"]["fault"]["fault"], partialsquad["trac"]["fault"]["fault"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)

    # Observation coordinates for plotting
    npts = 300
    obswidth = 5e3
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
    quadxslip = zeros(2, nels, 3)
    quadyslip = zeros(2, nels, 3)
    quadxslip[1, :, :] .= 1
    quadxslip[2, :, :] = transpose(reshape(slope .* xnodes, 3, nels))

    # Slip for partials
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
    for i in 2:size(slipconst)[1]
        # On fault evaluation
        dispconst = partialsconst["disp"]["fault"]["fault"] * slipconst[i, :]
        stressconst = partialsconst["stress"]["fault"]["fault"] * slipconst[i, :]
        tracconst = partialsconst["trac"]["fault"]["fault"] * slipconst[i, :]
        dispquad = partialsquad["disp"]["fault"]["fault"] * slipquad[i, :]
        stressquad = partialsquad["stress"]["fault"]["fault"] * slipquad[i, :]
        tracquad = partialsquad["trac"]["fault"]["fault"] * slipquad[i, :]

        # Far-field evaluation
        dispconstfarfield, stressconstfarfield = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], constxslip[i, :], constyslip[i, :], mu, nu)
        dispquadfarfield, stressquadfarfield = quaddispstress(slip2dispstress, xobs, yobs, els, idx["fault"], quadxslip[i, :, :], quadyslip[i, :, :], mu, nu)

        # Plot 9 panel results
        plotfunction(els, xcenters, ycenters, xnodes, ynodes, dispconst, stressconst, tracconst, dispquad, stressquad, tracquad, reshape(xobs, npts, npts), reshape(yobs, npts, npts), dispconstfarfield, stressconstfarfield, dispquadfarfield, stressquadfarfield)
    end

end
fig_constquadstresscomp()
