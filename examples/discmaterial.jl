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
    fontsize = 12
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
    discmaterial()

An attmept at the Crouch and Starfield level annulus solution
with varying material properties
"""
function discmaterial()
    close("all")
    mu = 3e10
    nu = 0.25
    p = -1.0e5 # Applied radial pressure over arc
    theta0 = deg2rad(1.0) # Arc length over which pressure is applied
    nels = 360
    R = nels / (2 * pi)
    npts = 50
    x, y = obsgrid(-R, -R, R, R, npts)
    r = @. sqrt(x^2 + y^2)

    # Define BEM Geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, nels)
    addelsez!(els, x1, y1, x2, y2, "a")
    idx = getidxdict(els)

    # Apply normal tractions everywhere and convert from radial to Cartesian
    xtracC = zeros(els.endidx)
    ytracC = zeros(els.endidx)
    thetaels = @. atan(els.ycenter[1:1:els.endidx], els.xcenter[1:1:els.endidx])
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtracC[i], ytracC[i] = els.rotmat[i, :, :] * normalTractions
    end

    #Kernels
    T_a_a, H_a_a = PUTC(slip2dispstress, els, idx["a"], idx["a"], mu, nu)
    Ueff = inv(H_a_a) * interleave(xtracC, ytracC)
    U, S = constdispstress(slip2dispstress, x, y, els, idx["a"], Ueff[1:2:end], Ueff[2:2:end], mu, nu)
    
    #! Summary figure
    figure(figsize=(30,20))
    fontsize = 12
    nrows = 3
    ncols = 3

    # Applied tractions
    subplot(nrows, ncols, 1)
    plot(rad2deg.(thetaels), xtracC, ".r")
    plot(rad2deg.(thetaels), ytracC, "+b")
    title("applied tractions", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    # Effective displacements
    subplot(nrows, ncols, 2)
    plot(rad2deg.(thetaels), Ueff[1:2:end], ".r")
    plot(rad2deg.(thetaels), Ueff[2:2:end], "+b")
    title("Effective displacements", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    # BEM solutions
    circle_subplot(nrows, ncols, 4, x, y, S[:, 1], npts, R, theta0, "Sxx (DDM)")
    circle_subplot(nrows, ncols, 5, x, y, S[:, 2], npts, R, theta0, "Syy (DDM)")
    circle_subplot(nrows, ncols, 6, x, y, S[:, 3], npts, R, theta0, "Sxy (DDM)")
end
discmaterial()
