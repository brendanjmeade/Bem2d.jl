using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d


"""
    discretized_arc(θstart, θend, radius, n_pts)

Generate regularly spaced eleemnts along an curved arc.
"""
function discretized_arc(θstart, θend, radius, n_pts)
    # Create geometry of discretized arc
    θrange = collect(LinRange(θstart, θend, n_pts + 1))
    x = @. radius * cos(θrange)
    y = @. radius * sin(θrange)
    x1 = x[1:1:end-1]
    x2 = x[2:1:end]
    y1 = y[1:1:end-1]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end


"""
    circle_subplot(nrows, ncols, plotidx, x, y, mat, npts, R, theta0, title_string)

Plot field (displacement, stress) within a circular disk and style
"""
function circle_subplot(nrows, ncols, plotidx, x, y, mat, npts, R, theta0, title_string)
    fontsize = 20
    contour_levels = 100
    contour_color = "white"
    contour_linewidth = 0.5
    color_scale = 1

    subplot(nrows, ncols, plotidx)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title(title_string, fontsize=fontsize)

    # Draw entire circle and region of applied tractions
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, 360)
    plot([x1, x2], [y1, y2], "-k", linewidth=2)
    x1, y1, x2, y2 = discretized_arc(-theta0, theta0, R, 50)
    plot([x1, x2], [y1, y2], "-r", linewidth=2)
    x1, y1, x2, y2 = discretized_arc(-theta0+deg2rad(180), theta0+deg2rad(180), R, 50)
    plot([x1, x2], [y1, y2], "-r", linewidth=2)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
end



"""
    discannulus()

Experimenting with material property variation with a disc embedded within 
an annulus
"""
function discannulus()
    close("all")
    mu = 3e10
    nu = 0.25
    p = -1.0e5 # Applied radial pressure over arc
    theta0 = deg2rad(45.0) # Arc length over which pressure is applied
    nels = 36
    r1 = 5 # radius of disc
    r2 = 10 # outer radius of annulus
    npts = 50
    x, y = Bem2d.obsgrid(-r2, -r2, r2, r2, npts)
    r = @. sqrt(x^2 + y^2)

    #! Define boundaries
    els = Bem2d.Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), r2, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "r2"
    end
    Bem2d.standardize_elements!(els)

    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), r1, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "r1"
    end
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    #! Apply normal displacement boundary conditions
    Ur2x = zeros(length(idx["r2"]))
    Ur2y = zeros(length(idx["r2"]))
    thetaels = @. atan(els.ycenter[idx["r2"]], els.xcenter[idx["r2"]])
    for i in 1:length(idx["r2"]) # Calcuate the x and y components of the displacements
        Ur2normal = [0; p] # Pressure in fault normal component only.
        Ur2x[idx["r2"][i]], Ur2y[idx["r2"][i]] = els.rotmat[idx["r2"][i], :, :] * Ur2normal
    end

    #! Zero out the tractions on the area without contact
    deleteidx = findall(x -> (x>theta0 && x<deg2rad(180)-theta0), thetaels)
    Ur2x[deleteidx] .= 0
    Ur2y[deleteidx] .= 0
    deleteidx = findall(x -> (x<-theta0 && x>-deg2rad(180)+theta0), thetaels)
    Ur2x[deleteidx] .= 0
    Ur2y[deleteidx] .= 0

    #! Kernels, T*: displacement to displacement, H*: displacement to traction
    Tstarr1, _, Hstarr1 = partialsconstdispstress(slip2dispstress, els, idx["r1"], idx["r1"], mu, nu)
    Tstarr2, _, Hstarr2 = partialsconstdispstress(slip2dispstress, els, idx["r2"], idx["r2"], mu, 0.5 * nu)
    Tr1Q, _, Hr1Q = partialsquaddispstress(slip2dispstress, els, idx["r1"], idx["r1"], mu, nu)
    Tr2Q, _, Hr2Q = partialsquaddispstress(slip2dispstress, els, idx["r2"], idx["r2"], mu, 0.5 * nu)


    @show cond(Tstarr1)
    @show cond(Tstarr2)
    @show cond(Tstarr1 - Tstarr2)
    # @infiltrate

    #! Internal stresses from applied displacements
    Udisp, Sdisp = constdispstress(slip2dispstress, x, y, els, idx["r2"], Ur2x, Ur2y, mu, nu)

    #! Isolate the values inside the circle
    nanidx = findall(x -> x > r2, r)
    Sdisp[nanidx, :] .= NaN
    
    #! Summary figure
    nrows = 3
    ncols = 3
    fontsize = 20
    figure(figsize=(20, 20))
    subplot(nrows, ncols, 1)
    # plot(rad2deg.(thetaels), Ur2x, ".r")
    # plot(rad2deg.(thetaels), Ur2y, "+b")
    title("applied displacements @ r2", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    # BEM solutions
    circle_subplot(nrows, ncols, 4, x, y, Sdisp[:, 1], npts, r2, theta0, "Sxx (DDM)")
    circle_subplot(nrows, ncols, 5, x, y, Sdisp[:, 2], npts, r2, theta0, "Syy (DDM)")
    circle_subplot(nrows, ncols, 6, x, y, Sdisp[:, 3], npts, r2, theta0, "Sxy (DDM)")

    tight_layout()
    show()
end
discannulus()
