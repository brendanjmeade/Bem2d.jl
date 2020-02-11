using Revise
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function local_subplot(x, y, mat, npts, title_string)
    fontsize = 20
    contour_levels = 50
    contour_color = "white"
    contour_linewidth = 0.5
    color_scale = 1e6
    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = PyPlot.colorbar(fraction=0.020, pad=0.05, extend = "both")
    PyPlot.contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    PyPlot.title(title_string, fontsize=fontsize)
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])
    PyPlot.xticks([])
    PyPlot.yticks([])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
end

function ex_flamant()
    mu = 3e10
    nu = 0.25
    npts = 200
    obswidth = 1000
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    r = @. sqrt(x^2 + y^2)
    θ = @. rad2deg(atan(y, x))
    fx = 0.0
    fy = 1.0
    σrr = @. -2.0/(pi*r) * (fx*cosd(θ) + fy*sind(θ))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))

    # Convert cylindrical stresses to Cartesian
    # Swap the transpose to go from Cartesian to cylindrical
    σxx_flamant = @. σrr * cosd(θ) * cosd(θ)
    σyy_flamant = @. σrr * sind(θ) * sind(θ)
    σxy_flamant = @. σrr * sind(θ) * cosd(θ)

    # TODO: Try the matrix rotation approach again
    # σxx = zeros(length(x))
    # σyy = zeros(length(x))
    # σxy = zeros(length(x))
    # for i in 1:length(x) # Project a single matrix to Cartesian coordinates
    #     cylindrical_stress_tensor = [σrr[i] σrθ[i] ; σrθ[i] σθθ[i]]
    #     transformation_matrix = [cosd(θ[i]) -sind(θ[i]) ; sind(θ[i]) cosd(θ[i])]
    #     cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transpose(transformation_matrix)
    #     σxx_flamant[i] = cartesian_stress_tensor[1, 1]
    #     σyy_flamant[i] = cartesian_stress_tensor[2, 2]
    #     σxy_flamant[i] = cartesian_stress_tensor[1, 2]
    # end



    # Try the Flamant solution from Crouch and Starfield section 3.1 (z-line load on "half-plane")
    σxx = @. -2*fy/pi * (x^2*y) / (x^2+y^2)^2
    σyy = @. -2*fy/pi * (y^3) / (x^2+y^2)^2
    σxy = @. -2*fy/pi * (x*y^2) / (x^2+y^2)^2

    # BEM solution
    els = Bem2d.Elements(Int(2))
    els.x1[1] = -0.5
    els.y1[1] = 0.0
    els.x2[1] = 0.5
    els.y2[1] = 0.0
    els.name[1] = "point"
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    # Given the traction boundary conditions calcuate the induced displacements on each element
    xdisp = zeros(els.endidx)
    ydisp = zeros(els.endidx)
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["point"][i], idx["point"][i], mu, nu)
        U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["point"][i], idx["point"][i], mu, nu)
        xdisp[i], ydisp[i] = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * [fx[i]; fy[i]]
    end
    dispall = Bem2d.interleave(xdisp, ydisp)

    # Streses from tractions
    _, stresstrac = Bem2d.constdispstress(trac2dispstress, x, y, els, idx["point"], fx, fy, mu, nu)

    # Stresses from traction induced displacements
    _, stressdisp = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["point"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    PyPlot.close("all")
    PyPlot.figure(figsize=(40,20))

    # Analytic Flamant
    PyPlot.subplot(3, 6, 1)
    local_subplot(x, y, σrr, npts, L"\sigma_{rr} \; \mathrm{(Wikipedia)}")
    PyPlot.subplot(3, 6, 2)
    local_subplot(x, y, σθθ, npts, L"\sigma_{\theta\theta} \; \mathrm{(Wikipedia)}")
    PyPlot.subplot(3, 6, 3)
    local_subplot(x, y, σrθ, npts, L"\sigma_{r\theta} \; \mathrm{(Wikipedia)}")
    PyPlot.subplot(3, 6, 7)
    local_subplot(x, y, σxx_flamant, npts, L"\sigma_{xx} \; \mathrm{(Wikipedia \; Cartesian)}")
    PyPlot.subplot(3, 6, 8)
    local_subplot(x, y, σyy_flamant, npts, L"\sigma_{yy} \; \mathrm{(Wikipedia \; Cartesian)}")
    PyPlot.subplot(3, 6, 9)
    local_subplot(x, y, σxy_flamant, npts, L"\sigma_{xy} \; \mathrm{(Wikipedia \; Cartesian)}")
    PyPlot.subplot(3, 6, 13)
    local_subplot(x, y, σxx, npts, L"\sigma_{xx} \; \mathrm{(CS \; Cartesian)}")
    PyPlot.subplot(3, 6, 14)
    local_subplot(x, y, σyy, npts, L"\sigma_{yy} \; \mathrm{(CS \; Cartesian)}")
    PyPlot.subplot(3, 6, 15)
    local_subplot(x, y, σxy, npts, L"\sigma_{xy} \; \mathrm{(CS \; Cartesian)}")

    # BEM Flamant
    PyPlot.subplot(3, 6, 4)
    local_subplot(x, y, stresstrac[:, 1], npts, L"\sigma_{xx} \; \mathrm{(traction)}")
    PyPlot.subplot(3, 6, 5)
    local_subplot(x, y, stresstrac[:, 2], npts, L"\sigma_{yy} \; \mathrm{(traction)}")
    PyPlot.subplot(3, 6, 6)
    local_subplot(x, y, stresstrac[:, 3], npts, L"\sigma_{xy} \; \mathrm{(traction)}")
    PyPlot.subplot(3, 6, 10)
    local_subplot(x, y, stressdisp[:, 1], npts, L"\sigma_{xx} \; \mathrm{(displacement)}")
    PyPlot.subplot(3, 6, 11)
    local_subplot(x, y, stressdisp[:, 2], npts, L"\sigma_{yy} \; \mathrm{(displacement)}")
    PyPlot.subplot(3, 6, 12)
    local_subplot(x, y, stressdisp[:, 3], npts, L"\sigma_{xy} \; \mathrm{(displacement)}")
    PyPlot.subplot(3, 6, 16)
    local_subplot(x, y, stresstrac[:, 1]+stressdisp[:, 1], npts, L"\sigma_{xx} \; \mathrm{(total)}")
    PyPlot.subplot(3, 6, 17)
    local_subplot(x, y, stresstrac[:, 2]+stressdisp[:, 2], npts, L"\sigma_{yy} \; \mathrm{(total)}")
    PyPlot.subplot(3, 6, 18)
    local_subplot(x, y, stresstrac[:, 3]+stressdisp[:, 3], npts, L"\sigma_{xy} \; \mathrm{(total)}")

    PyPlot.tight_layout()
    PyPlot.show()
end
ex_flamant()
