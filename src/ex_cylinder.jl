using Revise
using LaTeXStrings
using PyPlot
using Infiltrator
using Bem2d

function ex_cylinder()
    # Observation coordinates for far-field calculation
    npts = 200
    obswidth = 2e3
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Convert Cartesian to cylindrical coordinates
    r = @. sqrt(x^2 + y^2)
    θ = @. rad2deg(atan(y, x))

    # Solution from Hondros (1959) as summarized by Wei and Chau 2013
    σrr = zeros(length(x))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))
    p = 1.0e5 # Applied radial pressure over arc
    θ0 = 10.0 # Arc length over which pressure is applied
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
    # http://solidmechanics.org/text/AppendixD/AppendixD.htm
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

    # Plot contours with Makie/AbstractPlotting
    fontsize = 24
    contour_levels = 20
    contour_color = "black"
    contour_linewidth = 1.0
    PyPlot.close("all")
    PyPlot.figure(figsize=(30,15))

    PyPlot.subplot(2, 3, 1)
    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σrr, npts, npts), levels=contour_levels)
    PyPlot.colorbar()
    PyPlot.contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σrr, npts, npts), levels=contour_levels, colors="black")
    PyPlot.title(L"\sigma_{rr}", fontsize=fontsize)
    PyPlot.gca().set_aspect("equal")
    PyPlot.colorbar()

    PyPlot.subplot(2, 3, 2)
    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σθθ, npts, npts), levels=contour_levels)
    PyPlot.colorbar()
    PyPlot.contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σθθ, npts, npts), levels=contour_levels, colors="black")
    PyPlot.title(L"\sigma_{\theta\theta}", fontsize=fontsize)
    PyPlot.gca().set_aspect("equal")

    PyPlot.subplot(2, 3, 3)
    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σrθ, npts, npts), 20)
    PyPlot.title(L"\sigma_{r\theta}", fontsize=fontsize)
    PyPlot.gca().set_aspect("equal")
    PyPlot.colorbar()

    PyPlot.subplot(2, 3, 4)
    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σxx, npts, npts), 20)
    PyPlot.title(L"\sigma_{xx}", fontsize=fontsize)
    PyPlot.gca().set_aspect("equal")
    PyPlot.colorbar()

    PyPlot.subplot(2, 3, 5)
    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σyy, npts, npts), 20)
    PyPlot.title(L"\sigma_{yy}", fontsize=fontsize)
    PyPlot.gca().set_aspect("equal")
    PyPlot.colorbar()

    PyPlot.subplot(2, 3, 6)
    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(σxy, npts, npts), 20)
    PyPlot.title(L"\sigma_{xy}", fontsize=fontsize)
    PyPlot.gca().set_aspect("equal")
    PyPlot.colorbar()

    PyPlot.tight_layout()
    PyPlot.show()
end
ex_cylinder()
