using Revise
using LaTeXStrings
using PyPlot
using Infiltrator
using Bem2d

function circle_subplot(x, y, mat, npts, title_string)
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
    PyPlot.title(title_string, fontsize=fontsize)
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
    σrr = zeros(length(x))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))
    p = 1.0e5 # Applied radial pressure over arc
    θ0 = 70.0 # Arc length over which pressure is applied
    R = 1.0e3 # Radius of disc
    mmax = 10 # Max number of terms in Hondros series

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
    circle_subplot(x, y, σrr, npts, L"\sigma_{rr}")
    PyPlot.subplot(2, 3, 2)
    circle_subplot(x, y, σθθ, npts, L"\sigma_{\theta\theta}")
    PyPlot.subplot(2, 3, 3)
    circle_subplot(x, y, σrθ, npts, L"\sigma_{r\theta}")
    PyPlot.subplot(2, 3, 4)
    circle_subplot(x, y, σxx, npts, L"\sigma_{xx}")
    PyPlot.subplot(2, 3, 5)
    circle_subplot(x, y, σyy, npts, L"\sigma_{yy}")
    PyPlot.subplot(2, 3, 6)
    circle_subplot(x, y, σxy, npts, L"\sigma_{xy}")
    PyPlot.tight_layout()
    PyPlot.show()
end
ex_cylinder()
