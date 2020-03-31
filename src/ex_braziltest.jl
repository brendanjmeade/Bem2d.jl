using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

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

function circle_subplot(nrows, ncols, plotidx, x, y, mat, npts, R, θ0, title_string)
    fontsize = 20
    contour_levels = 20
    contour_color = "white"
    contour_linewidth = 0.5
    color_scale = 1

    subplot(nrows, ncols, plotidx)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.020, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    # cbar.set_label(label=title_string * " (Pa)", fontsize=fontsize)
    PyPlot.contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    PyPlot.title(title_string, fontsize=fontsize)
    # Draw entirre circle
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, 360)
    for i in 1:length(x1) # Draw arc where compression is being applied
        plot([x1[i], x2[i]], [y1[i], y2[i]], "-k", linewidth=2)
    end

    # Draw right-hand side of applied compression arc
    x1, y1, x2, y2 = discretized_arc(-θ0, θ0, R, 50)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-r", linewidth=2)
    end

    # Draw left-hand side of applied compression arc
    x1, y1, x2, y2 = discretized_arc(-θ0+deg2rad(180), θ0+deg2rad(180), R, 50)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-r", linewidth=2)
    end

    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
end

function normalizenan(vec)
    vec = vec .- minimum(filter(!isnan, vec))
    vec = vec ./ maximum(filter(!isnan, vec))
    return vec
end

function ex_braziltest()
    PyPlot.close("all")
    mu = 3e10
    nu = 0.25

    # Solution from Hondros (1959) as summarized by Wei and Chau 2013
    p = -1.0e5 # Applied radial pressure over arc
    θ0 = deg2rad(1.0) # Arc length over which pressure is applied
    nels = 360
    R = nels / (2 * pi) # Radius of disc ensuring that all elements are 1 unit long
    mmax = 1000 # Max number of terms in Hondros series
    npts = 50
    x, y = Bem2d.obsgrid(-R, -R, R, R, npts)
    r = @. sqrt(x^2 + y^2)
    θ = @. atan(y, x)

    #! Analytic stresses in cylindrical coordinates
    σrr = zeros(length(x))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))
    σrrconstterm = 2.0 * θ0 * -p / deg2rad(180)
    σθθconstterm = 2.0 * θ0 * -p / deg2rad(180)
    leadingterm = 2.0 * -p / deg2rad(180)
    for m in 1:mmax
        σrr += @. (r/R)^(2*m-2) * (1-(1-1/m)*(r/R)^2) * sin(2*m*θ0) * cos(2*m*θ)
        σθθ += @. (r/R)^(2*m-2) * (1-(1+1/m)*(r/R)^2) * sin(2*m*θ0) * cos(2*m*θ)
        σrθ += @. ((r/R)^(2*m) - (r/R)^(2*m-2)) * sin(2*m*θ0) * sin(2*m*θ)
    end
    σrr = @. σrrconstterm + leadingterm * σrr
    σθθ = @. σθθconstterm - leadingterm * σθθ
    σrθ = @. leadingterm * σrθ

    #! Convert analytic cylindrical stresses to Cartesian
    σxx = zeros(length(x))
    σyy = zeros(length(x))
    σxy = zeros(length(x))
    for i in 1:length(x) # Project a single matrix to Cartesian coordinates
        cylindrical_stress_tensor = [σrr[i] σrθ[i] ; σrθ[i] σθθ[i]]
        transformation_matrix = [cos(θ[i]) -sin(θ[i]) ; sin(θ[i]) cos(θ[i])]
        cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transpose(transformation_matrix)
        σxx[i] = cartesian_stress_tensor[1, 1]
        σyy[i] = cartesian_stress_tensor[2, 2]
        σxy[i] = cartesian_stress_tensor[1, 2]
    end

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
    partialsconst = Bem2d.initpartials(els)

    #! Apply normal tractions everywhere and convert from radial to Cartesian
    xtrac = zeros(els.endidx)
    ytrac = zeros(els.endidx)
    θels = @. atan(els.ycenter[1:1:els.endidx], els.xcenter[1:1:els.endidx])
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtrac[i], ytrac[i] = els.rotmat[i, :, :] * normalTractions
        ytrac[i] = -ytrac[i] # Just a desperate attempt
    end

    #! Zero out the tractions on the area without contact
    deleteidx = findall(x -> (x>θ0 && x<deg2rad(180)-θ0), θels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0
    deleteidx = findall(x -> (x<-θ0 && x>-deg2rad(180)+θ0), θels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0

    #! Scale tractions by element lengths
    xtracscaled = xtrac ./ els.length[1:1:els.endidx]
    ytracscaled = ytrac ./ els.length[1:1:els.endidx]

    #! Kernels
    # U*: traction to displacement
    # T*: displacement to displacement
    # A*: traction to traction
    # H*: displacement to traction
    Tstar, Sstar, Hstar = partialsconstdispstress(slip2dispstress, els, idx["circle"], idx["circle"], mu, nu)
    Ustar, Dstar, Astar = partialsconstdispstress(trac2dispstress, els, idx["circle"], idx["circle"], mu, nu)

    #! Induced displacements from Direct BEM
    dispall = (inv(Tstar + 0.5 * I(size(Tstar)[1]))) * Ustar * interleave(xtrac, ytrac)

    #! Displacement discontinuity method (indirect)
    dispall = Ustar * interleave(xtrac, ytrac)

    #! Streses from tractions
    _, Strac = constdispstress(trac2dispstress, x, y, els, idx["circle"], xtracscaled, ytracscaled, mu, nu)

    #! Stresses from traction induced displacements
    _, Sdisp = constdispstress(slip2dispstress, x, y, els, idx["circle"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    # #! Isolate the values inside the circle
    to_nan_idx = findall(x -> x > 0.9 * R, r)
    σrr[to_nan_idx] .= NaN
    σθθ[to_nan_idx] .= NaN
    σrθ[to_nan_idx] .= NaN
    σxx[to_nan_idx] .= NaN
    σyy[to_nan_idx] .= NaN
    σxy[to_nan_idx] .= NaN
    Strac[to_nan_idx, 1] .= NaN
    Strac[to_nan_idx, 2] .= NaN
    Strac[to_nan_idx, 3] .= NaN
    Sdisp[to_nan_idx, 1] .= NaN
    Sdisp[to_nan_idx, 2] .= NaN
    Sdisp[to_nan_idx, 3] .= NaN
    Sbem = @. Strac + Sdisp

    #! Summary figure
    fontsize = 24
    figure(figsize=(30,20))
    nrows = 3
    ncols = 6

    #! Traction forcing
    subplot(3, 6, 1)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k")
    end
    quiver(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], xtrac, ytrac)
    title("applied tractions", fontsize=fontsize)
    xticks([])
    yticks([])
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    #! Induced displacements
    subplot(3, 6, 2)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k")
    end
    quiver(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], dispall[1:2:end], dispall[2:2:end], color="green")
    title("induced displacements", fontsize=fontsize)
    xticks([])
    yticks([])
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    #! Analytic solutions
    circle_subplot(nrows, ncols, 7, x, y, σrr, npts, R, θ0, "Srr (analytic}")
    circle_subplot(nrows, ncols, 8, x, y, σθθ, npts, R, θ0, "Sthetatheta (analytic)")
    circle_subplot(nrows, ncols, 9, x, y, σrθ, npts, R, θ0, "Srtheta(analytic)")
    circle_subplot(nrows, ncols, 13, x, y, σxx, npts, R, θ0, "Sxx (analytic)")
    circle_subplot(nrows, ncols, 14, x, y, σyy, npts, R, θ0, "Syy (analytic)")
    circle_subplot(nrows, ncols, 15, x, y, σxy, npts, R, θ0, "Sxy (analytic)")

    #! BEM solutions
    circle_subplot(nrows, ncols, 4, x, y, Strac[:, 1], npts, R, θ0, "Sxx (tractions only)")
    circle_subplot(nrows, ncols, 5, x, y, Strac[:, 2], npts, R, θ0, "Syy (tractions only)")
    circle_subplot(nrows, ncols, 6, x, y, Strac[:, 3], npts, R, θ0, "Sxy (tractions only)")
    circle_subplot(nrows, ncols, 10, x, y, Sdisp[:, 1], npts, R, θ0, "Sxx (DDM)")
    circle_subplot(nrows, ncols, 11, x, y, Sdisp[:, 2], npts, R, θ0, "Syy (DDM)")
    circle_subplot(nrows, ncols, 12, x, y, Sdisp[:, 3], npts, R, θ0, "Sxy (DDM)")
    circle_subplot(nrows, ncols, 16, x, y, Sbem[:, 1], npts, R, θ0, "Sxx( sum)")
    circle_subplot(nrows, ncols, 17, x, y, Sbem[:, 2], npts, R, θ0, "Syy (sum)")
    circle_subplot(nrows, ncols, 18, x, y, Sbem[:, 3], npts, R, θ0, "Sxy (sum)")
    tight_layout()
    show()
end
ex_braziltest()
