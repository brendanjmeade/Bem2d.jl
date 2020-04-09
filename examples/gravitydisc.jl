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
    gravitydisc()

Experiments with gravity body force.
"""
function gravitydisc()
    close("all")
    mu = 3e10
    nu = 0.25
    p = -1.0e5 # Applied radial pressure over arc
    theta0 = deg2rad(15.0) # Arc length over which pressure is applied
    nels = 360
    R = nels / (2 * pi)
    npts = 50
    x, y = Bem2d.obsgrid(-R, -R, R, R, npts)
    r = @. sqrt(x^2 + y^2)

    #! BEM solution
    els = Bem2d.Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "circle"
    end
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    #! Apply normal tractions everywhere and convert from radial to Cartesian
    xtracC = zeros(els.endidx)
    ytracC = zeros(els.endidx)
    # xtracC[30] = 1.0
    thetaels = @. atan(els.ycenter[1:1:els.endidx], els.xcenter[1:1:els.endidx])
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtracC[i], ytracC[i] = els.rotmat[i, :, :] * normalTractions
    end

    #! Zero out the tractions on the area without contact
    deleteidx = findall(x -> (x>theta0 && x<deg2rad(180)-theta0), thetaels)
    xtracC[deleteidx] .= 0
    ytracC[deleteidx] .= 0
    deleteidx = findall(x -> (x<-theta0 && x>-deg2rad(180)+theta0), thetaels)
    xtracC[deleteidx] .= 0
    ytracC[deleteidx] .= 0

    #! Kernels, T*: displacement to displacement, H*: displacement to traction
    @time TstarC, _, HstarC = partialsconstdispstress(slip2dispstress, els, idx["circle"], idx["circle"], mu, nu)

    #! CONSTANT CASE
    #! Applied tractions -> effective displacements -> internal stresses
    DeffC = TstarC * interleave(xtracC, ytracC)
    # @time _, SdispC = constdispstress(slip2dispstress, x, y, els, idx["circle"], DeffC[1:2:end], DeffC[2:2:end], mu, nu)
    @time _, SdispC = constdispstress(slip2dispstress, x, y, els, idx["circle"], xtracC, ytracC, mu, nu)


    #! Isolate the values inside the circle
    nanidx = findall(x -> x > R, r)
    SdispC[nanidx, :] .= NaN
    # @infiltrate

    #! 9-panel plot for quadratic only
    figure(figsize=(15,15))
    nrows = 1
    ncols = 3

    # BEM solutions
    circle_subplot(nrows, ncols, 1, x, y, SdispC[:, 1], npts, R, theta0, L"\sigma_{xx} \; \mathrm{(BEM)}")
    circle_subplot(nrows, ncols, 2, x, y, SdispC[:, 2], npts, R, theta0, L"\sigma_{yy} \; \mathrm{(BEM)}")
    circle_subplot(nrows, ncols, 3, x, y, SdispC[:, 3], npts, R, theta0, L"\sigma_{xy} \; \mathrm{(BEM)}")
    tight_layout()
    show()
    @infiltrate
end
gravitydisc()
