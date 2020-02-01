using Revise
using LaTeXStrings
using PyPlot
using Infiltrator
using Bem2d

function discretized_arc(θstart, θend, radius, n_pts)
    # Create geometry of discretized arc
    θrange = collect(LinRange(θstart, θend, n_pts))
    x = @. radius * cosd(θrange)
    y = @. radius * sind(θrange)
    x1 = x[1:1:end-1]
    x2 = x[2:1:end]
    y1 = y[1:1:end-1]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end

function circle_subplot(x, y, mat, npts, R, θ0, title_string)
    fontsize = 20
    contour_levels = 10
    contour_color = "black"
    contour_linewidth = 0.5
    color_scale = 1e6
    mat = @. mat / 1e6 # Convert from Pascals to mega Pascals

    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = PyPlot.colorbar(fraction=0.020, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.set_label(label=title_string * " (MPa)", fontsize=fontsize)
    PyPlot.contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)

    x1, y1, x2, y2 = discretized_arc(-180, 180, R, 360)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-k", linewidth=10)
    end

    x1, y1, x2, y2 = discretized_arc(-θ0, θ0, R, 50)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-r", linewidth=5)
    end

    x1, y1, x2, y2 = discretized_arc(-θ0+180, θ0+180, R, 50)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-r", linewidth=5)
    end

    # PyPlot.title(title_string, fontsize=fontsize)
    PyPlot.xlabel(L"x (m)", fontsize=fontsize)
    PyPlot.ylabel(L"y (m)", fontsize=fontsize)
    PyPlot.xlim([-1100, 1100])
    PyPlot.ylim([-1100, 1100])
    PyPlot.xticks([-1000, 0, 1000])
    PyPlot.yticks([-1000, 0, 1000])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
end

function ex_cylinder()
    # Observation coordinates for far-field calculation
    npts = 200
    obswidth = 1e3
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Convert Cartesian to cylindrical coordinates
    r = @. sqrt(x^2 + y^2)
    θ = @. rad2deg(atan(y, x))

    # Solution from Hondros (1959) as summarized by Wei and Chau 2013
    p = 1.0e5 # Applied radial pressure over arc
    θ0 = 70.0 # Arc length over which pressure is applied
    R = 1.0e3 # Radius of disc
    mmax = 10 # Max number of terms in Hondros series

    σrr = zeros(length(x))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))
    σrrconstterm = 2.0 * θ0 * p / π
    σθθconstterm = 2.0 * θ0 * p / π
    leadingterm = 2.0 * p / π
    for m in 1:mmax
        σrr += @. (r/R)^(2*m-2) * (1-(1-1/m)*(r/R)^2) * sind(2*m*θ0) * cosd(2*m*θ)
        σθθ += @. (r/R)^(2*m-2) * (1-(1+1/m)*(r/R)^2) * sind(2*m*θ0) * cosd(2*m*θ)
        σrθ += @. ((r/R)^(2*m) - (r/R)^(2*m-2)) * sind(2*m*θ0) * sind(2*m*θ)
    end
    σrr = @. σrrconstterm + leadingterm * σrr
    σθθ = @. σθθconstterm + leadingterm * σθθ
    σrθ = @. leadingterm * σrθ

    # Convert cylindrical stresses to Cartesian
    # Swap the transpose to go from Cartesian to cylindrical
    σxx = zeros(length(x))
    σyy = zeros(length(x))
    σxy = zeros(length(x))
    for i in 1:length(x) # Project a single matrix to Cartesian coordinates
        cylindrical_stress_tensor = [σrr[i] σrθ[i] ; σrθ[i] σθθ[i]]
        transformation_matrix = [cosd(θ[i]) -sind(θ[i]) ; sind(θ[i]) cosd(θ[i])]
        cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transpose(transformation_matrix)
        σxx[i] = cartesian_stress_tensor[1, 1]
        σyy[i] = cartesian_stress_tensor[2, 2]
        σxy[i] = cartesian_stress_tensor[1, 2]
    end

    # Try setting a few values to NaN and see if we can isolate the circle
    to_nan_idx = findall(x -> x > R, r)
    σrr[to_nan_idx] .= NaN
    σθθ[to_nan_idx] .= NaN
    σrθ[to_nan_idx] .= NaN
    σxx[to_nan_idx] .= NaN
    σyy[to_nan_idx] .= NaN
    σxy[to_nan_idx] .= NaN

    # Plot contours with PyPlot
    fontsize = 24
    contour_levels = 10
    contour_color = "black"
    contour_linewidth = 0.5
    PyPlot.close("all")
    PyPlot.figure(figsize=(30,15))
    PyPlot.subplot(2, 3, 1)
    circle_subplot(x, y, σrr, npts, R, θ0, L"\sigma_{rr}")
    PyPlot.subplot(2, 3, 2)
    circle_subplot(x, y, σθθ, npts, R, θ0, L"\sigma_{\theta\theta}")
    PyPlot.subplot(2, 3, 3)
    circle_subplot(x, y, σrθ, npts, R, θ0, L"\sigma_{r\theta}")
    PyPlot.subplot(2, 3, 4)
    circle_subplot(x, y, σxx, npts, R, θ0, L"\sigma_{xx}")
    PyPlot.subplot(2, 3, 5)
    circle_subplot(x, y, σyy, npts, R, θ0, L"\sigma_{yy}")
    PyPlot.subplot(2, 3, 6)
    circle_subplot(x, y, σxy, npts, R, θ0, L"\sigma_{xy}")
    PyPlot.suptitle(string(θ0), fontsize=30)
    PyPlot.tight_layout()
    PyPlot.show()
end
ex_cylinder()
