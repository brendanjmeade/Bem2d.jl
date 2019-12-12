using Revise
using PyCall
using PyPlot
using Bem2d

function plotfunction(els, xcenters, ycenters, xnodes, ynodes, uconst, σconst, tconst, uquad, σquad, tquad, xobs, yobs, uconstfarfield, σconstfarfield, uquadfarfield, σquadfarfield)
    scalestress = 1e6
    ncontours = 30
    fontsize = 5
    xmin = -20.0e3
    xmax = 20.0e3
    cmap = rycroftcmap()

    figure(figsize = (6, 5), dpi=300)
    ax = subplot(3, 3, 1)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{xx}$ CS (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape(σconstfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize = fontsize); xlim([xmin, xmax]); xticks([])
    yticks([-20000, 0, 20000]);
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 2)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 2] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{yy}$ CS (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape(σconstfarfield[:, 2] ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 3)
    contourf(xobs, yobs, reshape(σconstfarfield[:, 3] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{xy}$ CS (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape(σconstfarfield[:, 3] ./ scalestress, size(xobs)), ncontours,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 4)
    contourf(xobs, yobs, reshape(σquadfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{xx}$ 3NQ (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape(σquadfarfield[:, 1] ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize = fontsize); xlim([xmin, xmax]); xticks([])
    yticks([-20000, 0, 20000]);
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 5)
    contourf(xobs, yobs, reshape(σquadfarfield[:, 2] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{xy}$ 3NQ (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape(σquadfarfield[:, 2] ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 6)
    contourf(xobs, yobs, reshape(σquadfarfield[:, 3] ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{yy}$ 3NQ (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape(σquadfarfield[:, 3] ./ scalestress, size(xobs)), ncontours,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 7)
    contourf(xobs, yobs, reshape((σconstfarfield[:, 1] - σquadfarfield[:, 1]) ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{xx}$ 3NQ (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape((σconstfarfield[:, 1] - σquadfarfield[:, 1]) ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize = fontsize); xlim([xmin, xmax]); xticks([])
    yticks([-20000, 0, 20000]);
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 8)
    contourf(xobs, yobs, reshape((σconstfarfield[:, 2] - σquadfarfield[:, 2]) ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{xy}$ 3NQ (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape((σconstfarfield[:, 2] - σquadfarfield[:, 2]) ./ scalestress, size(xobs)), ncontours, linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 3, 9)
    contourf(xobs, yobs, reshape((σconstfarfield[:, 3] - σquadfarfield[:, 3]) ./ scalestress, size(xobs)), ncontours, cmap = cmap)
    cbar = colorbar(orientation = "horizontal", fraction = 0.20, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.set_label(label = L"$\sigma_{yy}$ 3NQ (MPa)", fontsize = fontsize)
    contour(xobs, yobs, reshape((σconstfarfield[:, 3] - σquadfarfield[:, 3]) ./ scalestress, size(xobs)), ncontours,  linewidths = 0.25, colors = "k", linestyles = "-")
    gca().set_aspect("equal")
    xlim([xmin, xmax]); xticks([]); yticks([])
    ax.tick_params(labelsize = fontsize)

    tight_layout()
    show()
end

function plotprofiles(xobs, yobs, uconstfarfield, σconstfarfield, uquadfarfield, σquadfarfield)
    fontsize = 10
    linewidth = 1.0
    xmin = -3.0e3
    xmax = 3.0e3
    ymin = 3
    ymax = 7
    
    figure(figsize = (5, 5), dpi=300)
    ax = subplot(3, 1, 1)
    plot(yobs, log10.(abs.(σconstfarfield[:, 1])), "-b", linewidth=linewidth)
    plot(yobs, log10.(abs.(σquadfarfield[:, 1])), "--r", linewidth=linewidth)
    xlabel(L"$y$ (m)", fontsize = fontsize);
    xlim([xmin, xmax])
    ylabel(L"$\sigma_{xx}$ (m)", fontsize = fontsize);
    xticks([xmin, 0, xmax]);
    ylim([ymin, ymax])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 1, 2)
    plot(yobs, log10.(abs.(σconstfarfield[:, 2])), "-b", linewidth=linewidth)
    plot(yobs, log10.(abs.(σquadfarfield[:, 2])), "--r", linewidth=linewidth)
    xlabel(L"$y$ (m)", fontsize = fontsize);
    xlim([xmin, xmax])
    ylabel(L"$\sigma_{yy}$ (m)", fontsize = fontsize);
    xticks([xmin, 0, xmax]);
    ylim([ymin, ymax])
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 1, 3)
    plot(yobs, log10.(abs.(σconstfarfield[:, 3])), "-b", linewidth=linewidth)
    plot(yobs, log10.(abs.(σquadfarfield[:, 3])), "--r", linewidth=linewidth)    
    xlabel(L"$y$ (m)", fontsize = fontsize);
    xlim([xmin, xmax])
    ylabel(L"$\log\|\sigma_{xy}\|$ (m)", fontsize = fontsize);
    xticks([xmin, 0, xmax]);
    ylim([ymin, ymax])
    ax.tick_params(labelsize = fontsize)

    tight_layout()
    show()
end

function fig_constquadstresscomp()
    # Material and geometric constants
    mu = 3e10
    nu = 0.25
    nels = 19
    els = Elements(Int(1e5))
    L = 10000
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nels)
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
    npts = 100
    obswidth = 20e3
    xobs, yobs = obsgrid(-obswidth, -20.0e3, obswidth, 20.0e3, npts)

    # Evaluation points
    xcenters = els.xcenter[1:els.endidx]
    ycenters = els.ycenter[1:els.endidx]
    xnodes = sort(els.xnodes[1:els.endidx, :][:])
    ynodes = sort(els.ynodes[1:els.endidx, :][:])

    # Slip for volumetric models
    slope = 0.0001
    constxslip = zeros(2, nels)
    constyslip = zeros(2, nels)
    constxslip[1, :] .= 5.0
    constxslip[2, :] = slope .* xcenters
    quadxslip = zeros(2, nels, 3)
    quadyslip = zeros(2, nels, 3)
    quadxslip[1, :, :] .= 5.0
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
    nprofile = 1000
    xobs = 100 .* ones(nprofile)
    yobs = collect(LinRange(-3e3, 3e3, nprofile))
    close("all")
    for i in 1:size(slipconst)[1]
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

        # Plot profiles
        # plotprofiles(xobs, yobs, dispconstfarfield, stressconstfarfield, dispquadfarfield, stressquadfarfield)
        
        # Plot 9 panel results
        # plotfunction(els, xcenters, ycenters, xnodes, ynodes, dispconst, stressconst, tracconst, dispquad, stressquad, tracquad, reshape(xobs, npts, npts), reshape(yobs, npts, npts), dispconstfarfield, stressconstfarfield, dispquadfarfield, stressquadfarfield)
    end

    i = 1
    dispconst = partialsconst["disp"]["fault"]["fault"] * slipconst[i, :]
    stressconst = partialsconst["stress"]["fault"]["fault"] * slipconst[i, :]
    tracconst = partialsconst["trac"]["fault"]["fault"] * slipconst[i, :]
    dispquad = partialsquad["disp"]["fault"]["fault"] * slipquad[i, :]
    stressquad = partialsquad["stress"]["fault"]["fault"] * slipquad[i, :]
    tracquad = partialsquad["trac"]["fault"]["fault"] * slipquad[i, :]
    
    # Far-field evaluation
    dispconstfarfield, stressconstfarfield = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], constxslip[i, :], constyslip[i, :], mu, nu)
    dispquadfarfield, stressquadfarfield = quaddispstress(slip2dispstress, xobs, yobs, els, idx["fault"], quadxslip[i, :, :], quadyslip[i, :, :], mu, nu)
    
    fontsize = 10
    linewidth = 2.0
    xmin = -3.0e3
    xmax = 3.0e3
    ymin = 2
    ymax = 7
    
    figure(figsize = (7, 6), dpi=300)
    ax = subplot(3, 2, 1)
    plot(yobs, log10.(abs.(stressconstfarfield[:, 1])), "-b", linewidth=linewidth, label="CS")
    plot(yobs, log10.(abs.(stressquadfarfield[:, 1])), "--r", linewidth=linewidth, label="3QN")
    # xlabel(L"$y$ (m)", fontsize = fontsize);
    xlim([xmin, xmax])
    ylabel(L"$\log\|\sigma_{xx}\|$ (Pa)", fontsize = fontsize);
    xticks([])
    ylim([ymin, ymax])
    legend(frameon=false, loc="lower left")
    title("constant slip")
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 2, 3)
    plot(yobs, log10.(abs.(stressconstfarfield[:, 2])), "-b", linewidth=linewidth, label="CS")
    plot(yobs, log10.(abs.(stressquadfarfield[:, 2])), "--r", linewidth=linewidth, label="3QN")
    # xlabel(L"$y$ (m)", fontsize = fontsize)
    xlim([xmin, xmax])
    ylabel(L"$\log\|\sigma_{yy}\|$ (Pa)", fontsize = fontsize)
    xticks([])
    ylim([ymin, ymax])
    legend(frameon=false, loc="lower left")
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 2, 5)
    plot(yobs, log10.(abs.(stressconstfarfield[:, 3])), "-b", linewidth=linewidth, label="CS")
    plot(yobs, log10.(abs.(stressquadfarfield[:, 3])), "--r", linewidth=linewidth, label="3QN")
    xlabel(L"$y$ (m)", fontsize = fontsize)
    xlim([xmin, xmax])
    ylabel(L"$\log\|\sigma_{xy}\|$ (Pa)", fontsize = fontsize)
    xticks([xmin, 0, xmax])
    ylim([ymin, ymax])
    legend(frameon=false, loc="lower left")
    ax.tick_params(labelsize = fontsize)

    i = 2
    dispconst = partialsconst["disp"]["fault"]["fault"] * slipconst[i, :]
    stressconst = partialsconst["stress"]["fault"]["fault"] * slipconst[i, :]
    tracconst = partialsconst["trac"]["fault"]["fault"] * slipconst[i, :]
    dispquad = partialsquad["disp"]["fault"]["fault"] * slipquad[i, :]
    stressquad = partialsquad["stress"]["fault"]["fault"] * slipquad[i, :]
    tracquad = partialsquad["trac"]["fault"]["fault"] * slipquad[i, :]
    
    # Far-field evaluation
    dispconstfarfield, stressconstfarfield = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], constxslip[i, :], constyslip[i, :], mu, nu)
    dispquadfarfield, stressquadfarfield = quaddispstress(slip2dispstress, xobs, yobs, els, idx["fault"], quadxslip[i, :, :], quadyslip[i, :, :], mu, nu)
    
    ax = subplot(3, 2, 2)
    plot(yobs, log10.(abs.(stressconstfarfield[:, 1])), "-b", linewidth=linewidth, label="CS")
    plot(yobs, log10.(abs.(stressquadfarfield[:, 1])), "--r", linewidth=linewidth, label="3QN")
    # xlabel(L"$y$ (m)", fontsize = fontsize)
    xlim([xmin, xmax])
    # ylabel(L"$\log\|\sigma_{xx}\|$ (Pa)", fontsize = fontsize)
    xticks([])
    yticks([])
    ylim([ymin, ymax])
    title("linear slip")
    legend(frameon=false, loc="lower left")
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 2, 4)
    plot(yobs, log10.(abs.(stressconstfarfield[:, 2])), "-b", linewidth=linewidth, label="CS")
    plot(yobs, log10.(abs.(stressquadfarfield[:, 2])), "--r", linewidth=linewidth, label="3QN")
    # xlabel(L"$y$ (m)", fontsize = fontsize)
    xlim([xmin, xmax])
    # ylabel(L"$\log\|\sigma_{yy}\|$ (Pa)", fontsize = fontsize)
    xticks([])
    yticks([])
    ylim([ymin, ymax])
    legend(frameon=false, loc="lower left")
    ax.tick_params(labelsize = fontsize)

    ax = subplot(3, 2, 6)
    plot(yobs, log10.(abs.(stressconstfarfield[:, 3])), "-b", linewidth=linewidth, label="CS")
    plot(yobs, log10.(abs.(stressquadfarfield[:, 3])), "--r", linewidth=linewidth, label="3QN")
    xlabel(L"$y$ (m)", fontsize = fontsize);
    xlim([xmin, xmax])
    # ylabel(L"$\log\|\sigma_{xy}\|$ (Pa)", fontsize = fontsize)
    yticks([])
    xticks([xmin, 0, xmax])
    ylim([ymin, ymax])
    legend(frameon=false, loc="lower left")
    ax.tick_params(labelsize = fontsize)

    tight_layout()
    show()

    
end
fig_constquadstresscomp()
