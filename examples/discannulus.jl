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
    @time TstarQ, _, HstarQ = partialsquaddispstress(slip2dispstress, els, idx["circle"], idx["circle"], mu, nu)

    #! CONSTANT CASE
    #! Applied tractions -> effective displacements -> internal stresses
    DeffC = TstarC * interleave(xtracC, ytracC)
    # @time _, SdispC = constdispstress(slip2dispstress, x, y, els, idx["circle"], DeffC[1:2:end], DeffC[2:2:end], mu, nu)
    @time _, SdispC = constdispstress(slip2dispstress, x, y, els, idx["circle"], xtracC, ytracC, mu, nu)


    #! QUADRATIC CASE
    #! Applied tractions -> effective displacements -> internal stresses
    xtracQ = zeros(3 * length(xtracC))
    ytracQ = zeros(3 * length(ytracC))
    xtracQ[1:3:end] = xtracC
    xtracQ[2:3:end] = xtracC
    xtracQ[3:3:end] = xtracC
    ytracQ[1:3:end] = ytracC
    ytracQ[2:3:end] = ytracC
    ytracQ[3:3:end] = ytracC
    DeffQ = TstarQ * interleave(xtracQ, ytracQ)
    # @time _, SdispQ = quaddispstress(slip2dispstress, x, y, els, idx["circle"], quadstack(DeffQ[1:2:end]), quadstack(DeffQ[2:2:end]), mu, nu)
    @time _, SdispQ = quaddispstress(slip2dispstress, x, y, els, idx["circle"], quadstack(DeffQ[1:2:end]), quadstack(DeffQ[2:2:end]), mu, nu)


    #! Isolate the values inside the circle
    nanidx = findall(x -> x > R, r)
    SdispC[nanidx, :] .= NaN
    SdispQ[nanidx, :] .= NaN
    SresidualC = SdispC 
    SresidualQ = SdispQ
    
    #! Summary figure
    figure(figsize=(30,20))
    fontsize = 20
    nrows = 4
    ncols = 6

    #! Constant dislocation case
    # Applied tractions
    subplot(nrows, ncols, 1)
    plot(rad2deg.(thetaels), xtracC, ".r")
    plot(rad2deg.(thetaels), ytracC, "+b")
    title("applied tractions", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    # Effective displacements
    subplot(nrows, ncols, 2)
    plot(rad2deg.(thetaels), DeffC[1:2:end], ".r")
    plot(rad2deg.(thetaels), DeffC[2:2:end], "+b")
    title("Effective displacements", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    # BEM solutions
    circle_subplot(nrows, ncols, 13, x, y, SdispC[:, 1], npts, R, theta0, "Sxx (DDM)")
    circle_subplot(nrows, ncols, 14, x, y, SdispC[:, 2], npts, R, theta0, "Syy (DDM)")
    circle_subplot(nrows, ncols, 15, x, y, SdispC[:, 3], npts, R, theta0, "Sxy (DDM)")
    circle_subplot(nrows, ncols, 19, x, y, SresidualC[:, 1], npts, R, theta0, "Sxx (residual)")
    circle_subplot(nrows, ncols, 20, x, y, SresidualC[:, 2], npts, R, theta0, "Syy (residual)")
    circle_subplot(nrows, ncols, 21, x, y, SresidualC[:, 3], npts, R, theta0, "Sxy (residual)")

    #! Quadratic case
    # Applied tractions
    subplot(nrows, ncols, 4)
    plot(xtracQ, ".r")
    plot(ytracQ, "+b")
    title("applied tractions", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    # Effective displacements
    subplot(nrows, ncols, 5)
    plot(DeffQ[1:2:end], ".r")
    plot(DeffQ[2:2:end], "+b")
    title("Effective displacements", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    # BEM solutions
    circle_subplot(nrows, ncols, 16, x, y, SdispQ[:, 1], npts, R, theta0, "Sxx (DDM)")
    circle_subplot(nrows, ncols, 17, x, y, SdispQ[:, 2], npts, R, theta0, "Syy (DDM)")
    circle_subplot(nrows, ncols, 18, x, y, SdispQ[:, 3], npts, R, theta0, "Sxy (DDM)")
    circle_subplot(nrows, ncols, 22, x, y, SresidualQ[:, 1], npts, R, theta0, "Sxx (residual)")
    circle_subplot(nrows, ncols, 23, x, y, SresidualQ[:, 2], npts, R, theta0, "Syy (residual)")
    circle_subplot(nrows, ncols, 24, x, y, SresidualQ[:, 3], npts, R, theta0, "Sxy (residual)")
    tight_layout()
    show()
end
discannulus()
