using Revise
using PyPlot
using Infiltrator
using Bem2d

function local_subplot(x, y, mat, npts, title_string)
    fontsize = 20
    contour_levels = 20
    contour_color = "black"
    contour_linewidth = 0.5
    color_scale = 1e6
    mat = @. mat / 1e6 # Convert from Pascals to mega Pascals

    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = PyPlot.colorbar(fraction=0.020, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.set_label(label=title_string * " (MPa)", fontsize=fontsize)
    PyPlot.contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    PyPlot.xlabel(L"x \; (m)", fontsize=fontsize)
    PyPlot.ylabel(L"y \; (m)", fontsize=fontsize)
    PyPlot.xlim([-1100, 1100])
    PyPlot.ylim([-1100, 1100])
    PyPlot.xticks([-1000, 0, 1000])
    PyPlot.yticks([-1000, 0, 1000])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
end

function ex_flamant()
    # Observation coordinates
    npts = 50
    obswidth = 1000
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    r = @. sqrt(x^2 + y^2)
    θ = @. rad2deg(atan(y, x))

    fx = 1.0
    fy = 0.0
    σrr = zeros(length(x))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))
    σrr = @. -2.0/(pi*r) * (fx*cosd(θ) + fy*sind(θ))

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

    PyPlot.figure(figsize=(30,20))
    PyPlot.subplot(3, 3, 1)
    local_subplot(x, y, σrr, npts, L"\sigma_{rr}")
    PyPlot.subplot(3, 3, 2)
    local_subplot(x, y, σθθ, npts, L"\sigma_{\theta\theta}")
    PyPlot.subplot(3, 3, 3)
    local_subplot(x, y, σrθ, npts, L"\sigma_{r\theta}")
    PyPlot.subplot(3, 3, 4)
    local_subplot(x, y, σxx, npts, L"\sigma_{xx}")
    PyPlot.subplot(3, 3, 5)
    local_subplot(x, y, σyy, npts, L"\sigma_{yy}")
    PyPlot.subplot(3, 3, 6)
    local_subplot(x, y, σxy, npts, L"\sigma_{xy}")

end
ex_flamant()
