using Revise
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function local_subplot(x, y, mat, npts, title_string)
    fontsize = 20
    contour_levels = 50
    contour_levels = collect(LinRange(-0.01, 0.01, 31))

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


function stylesubplots_local()
    gca().set_aspect("equal")
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])

    # gca().set_xticks([])
    # gca().set_yticks([])
    return nothing
end

function plot18_local(els, x, y, disp1, stress1, string1, disp2, stress2, string2, title_string)    
    # Set contour levels for displacements and stresses
    contourvecdispy = collect(LinRange(-5e-11, 5e-11, 51))
    contourvecdispx = collect(LinRange(-1e-12, 1e-12, 51))
    contourvecstress = collect(LinRange(-1e-2, 1e-2, 51))
    cmap = PyPlot.get_cmap("seismic")
    fontsize = 30
    PyPlot.figure(figsize=(30, 20))

    # Fields from first model
    PyPlot.subplot(3, 6, 1)
    PyPlot.quiver(x[:], y[:], disp1[:, 1], disp1[:, 2], units="width", color="b")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("u", fontsize=fontsize)

    PyPlot.subplot(3, 6, 2)
    PyPlot.contourf(x, y, reshape(disp1[:, 1], size(x)), contourvecdispx, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 1], size(x)), contourvecdispx, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("ux (" * string1 * ")", fontsize=fontsize)
    
    PyPlot.subplot(3, 6, 3)
    PyPlot.contourf(x, y, reshape(disp1[:, 2], size(x)), contourvecdispy, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 2], size(x)), contourvecdispy, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("uy (" * string1 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 4)
    PyPlot.contourf(x, y, reshape(stress1[:, 1], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress1[:, 1], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("sxx (" * string1 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 5)
    PyPlot.contourf(x, y, reshape(stress1[:, 2], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress1[:, 2], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("syy (" * string1 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 6)
    PyPlot.contourf(x, y, reshape(stress1[:, 3], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress1[:, 3], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("sxy (" * string1 * ")", fontsize=fontsize)

    # Fields from second model
    PyPlot.subplot(3, 6, 7)
    PyPlot.quiver(x[:], y[:], disp2[:, 1], disp2[:, 2], units="width", color="b")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("u", fontsize=fontsize)

    PyPlot.subplot(3, 6, 8)
    PyPlot.contourf(x, y, reshape(disp2[:, 1], size(x)), contourvecdispx, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp2[:, 1], size(x)), contourvecdispx, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("ux (" * string2 * ")", fontsize=fontsize)
    
    PyPlot.subplot(3, 6, 9)
    PyPlot.contourf(x, y, reshape(disp2[:, 2], size(x)), contourvecdispy, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp2[:, 2], size(x)), contourvecdispy, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("uy (" * string2 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 10)
    PyPlot.contourf(x, y, reshape(stress2[:, 1], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress2[:, 1], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("sxx (" * string2 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 11)
    PyPlot.contourf(x, y, reshape(stress2[:, 2], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress2[:, 2], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("syy (" * string2 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 12)
    PyPlot.contourf(x, y, reshape(stress2[:, 3], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress2[:, 3], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("sxy (" * string2 * ")", fontsize=fontsize)

    # Residuals
    PyPlot.subplot(3, 6, 13)
    PyPlot.quiver(x[:], y[:], disp1[:, 1]-disp2[:, 1], disp1[:, 2]-disp2[:, 2], units="width", color="b")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("u", fontsize=fontsize)

    PyPlot.subplot(3, 6, 14)
    PyPlot.contourf(x, y, reshape(disp1[:, 1]-disp2[:, 1], size(x)), contourvecdispx, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 1]-disp2[:, 1], size(x)), contourvecdispx, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("uy (" * string1 * "-" * string2 * ")", fontsize=fontsize)
    
    PyPlot.subplot(3, 6, 15)
    PyPlot.contourf(x, y, reshape(disp1[:, 2]-disp2[:, 2], size(x)), contourvecdispy, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 2]-disp2[:, 2], size(x)), contourvecdispy, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("uy (" * string1 * "-" * string2 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 16)
    PyPlot.contourf(x, y, reshape(stress1[:, 1]-stress2[:, 1], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress1[:, 1]-stress2[:, 1], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("sxx (" * string1 * "-" * string2 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 17)
    PyPlot.contourf(x, y, reshape(stress1[:, 2]-stress2[:, 2], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress1[:, 2]-stress2[:, 2], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("syy (" * string1 * "-" * string2 * ")", fontsize=fontsize)

    PyPlot.subplot(3, 6, 18)
    PyPlot.contourf(x, y, reshape(stress1[:, 3]-stress2[:, 3], size(x)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(stress1[:, 3]-stress2[:, 3], size(x)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("sxy (" * string1 * "-" * string2 * ")", fontsize=fontsize)

    PyPlot.suptitle(title_string, fontsize=fontsize)
    PyPlot.show()
end


function ex_flamant()
    obswidth = 1000

    #! Try wierd symmetric quadratic spacing with small elements i the middle
    nquad = 100
    x1 = zeros(3 + 2 * (nquad-1))
    x2 = zeros(3 + 2 * (nquad-1))
    quadspacing = collect(LinRange(1, sqrt(3*obswidth), nquad).^2)
    x1[1:1:nquad-1] = -reverse(quadspacing[2:1:end], dims=1)
    x1[nquad:1:nquad+2] = [-1.0 -0.5 0.5]
    x1[nquad+3:1:end] = quadspacing[1:1:end-1]
    x2[1:1:nquad-1] = -reverse(quadspacing[1:1:end-1], dims=1)
    x2[nquad:1:nquad+2] = [-0.5 0.5 1.0]
    x2[nquad+3:1:end] = quadspacing[2:1:end]
    y1 = zeros(size(x1))
    y2 = zeros(size(x1))
    mididx = Int(floor(length(x1)/2)+1)

    mu = 3e10
    nu = 0.25
    npts = 50
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    r = @. sqrt(x^2 + y^2)
    θ = @. rad2deg(atan(y, x))
    nels = length(x1)
    fx = zeros(nels)
    fy = zeros(nels)
    fy[mididx] = 1.0

    #! Analytic solution from Wikipedia in cylindrical coordinates
    σrr = @. -2.0/(pi*r) * (fx[mididx]*cosd(θ) + fy[mididx]*sind(θ))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))

    σxx_flamant = zeros(length(x))
    σyy_flamant = zeros(length(x))
    σxy_flamant = zeros(length(x))
    for i in 1:length(x) # Project a single matrix to Cartesian coordinates
        cylindrical_stress_tensor = [σrr[i] σrθ[i] ; σrθ[i] σθθ[i]]
        transformation_matrix = [cosd(θ[i]) -sind(θ[i]) ; sind(θ[i]) cosd(θ[i])]
        cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transpose(transformation_matrix)
        σxx_flamant[i] = cartesian_stress_tensor[1, 1]
        σyy_flamant[i] = cartesian_stress_tensor[2, 2]
        σxy_flamant[i] = cartesian_stress_tensor[1, 2]
    end

    #! BEM solution
    els = Bem2d.Elements(Int(1e5))
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "line"
    end
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    #! Given the traction boundary conditions calcuate the induced displacements on each element
    xdisp = zeros(els.endidx)
    ydisp = zeros(els.endidx)
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["line"][i], idx["line"][i], mu, nu)
        U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["line"][i], idx["line"][i], mu, nu)
        xdisp[i], ydisp[i] = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * [fx[i]; fy[i]]
    end
    displocal = Bem2d.interleave(xdisp, ydisp)

    # Given the other tractions calcuate the induced displacements on the boudaries
    T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["line"], idx["line"], mu, nu)
    U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["line"], idx["line"], mu, nu)
    dispall = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * Bem2d.interleave(fx, fy)

    @show cond(T + 0.5 * LinearAlgebra.I(size(T)[1]))
    # U[diagind(U)] .= 0
    # dispall = U * Bem2d.interleave(fx, fy) + displocal

    # Streses from tractions
    _, stresstrac = Bem2d.constdispstress(trac2dispstress, x, y, els, idx["line"], fx, fy, mu, nu)

    # Stresses from traction induced displacements
    _, stressdisp = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["line"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    # Displacement solutions (https://en.wikipedia.org/wiki/Flamant_solution)
    kappaplainstrain = 3-4*nu
    kappaplanestress = (3-nu)/(1+nu)
    
    #! Set everything in the upper half-plane to zero because that is the exterior solution
    # to_nan_idx = findall(x -> x > 0.0, y)
    # σrr[to_nan_idx] .= NaN
    # σθθ[to_nan_idx] .= NaN
    # σrθ[to_nan_idx] .= NaN
    # σxx_flamant[to_nan_idx] .= NaN
    # σyy_flamant[to_nan_idx] .= NaN
    # σxy_flamant[to_nan_idx] .= NaN
    # stresstrac[to_nan_idx, 1] .= NaN
    # stresstrac[to_nan_idx, 2] .= NaN
    # stresstrac[to_nan_idx, 3] .= NaN
    # stressdisp[to_nan_idx, 1] .= NaN
    # stressdisp[to_nan_idx, 2] .= NaN
    # stressdisp[to_nan_idx, 3] .= NaN

    #! 18 panel plot
    PyPlot.close("all")
    dispbem = zeros(length(x), 2)
    stressbem = @. stresstrac + stressdisp
    dispanalytic = zeros(length(x), 2)
    stressanalytic = [σxx_flamant σyy_flamant σxy_flamant] 
    xobs = reshape(x, npts, npts)
    yobs = reshape(y, npts, npts)
    plot18_local(els, xobs, yobs, dispbem, stressbem, "BEM", dispanalytic, stressanalytic, "analytic", "Flamant (BEM vs. analytic)")

    #! Hacky plotting
    PyPlot.figure(figsize=(40,20))
    #! Analytic Flamant
    PyPlot.subplot(3, 6, 13)
    local_subplot(x, y, σxx_flamant, npts, L"\sigma_{xx} \; \mathrm{(Wikipedia \; Cartesian)}")
    PyPlot.subplot(3, 6, 14)
    local_subplot(x, y, σyy_flamant, npts, L"\sigma_{yy} \; \mathrm{(Wikipedia \; Cartesian)}")
    PyPlot.subplot(3, 6, 15)
    local_subplot(x, y, σxy_flamant, npts, L"\sigma_{xy} \; \mathrm{(Wikipedia \; Cartesian)}")

    #! BEM Flamant
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
    local_subplot(x, y, stresstrac[:, 1] .- stressdisp[:, 1], npts, L"\sigma_{xx} \; \mathrm{(total)}")
    PyPlot.subplot(3, 6, 17)
    local_subplot(x, y, stresstrac[:, 2] .- stressdisp[:, 2], npts, L"\sigma_{yy} \; \mathrm{(total)}")
    PyPlot.subplot(3, 6, 18)
    local_subplot(x, y, stresstrac[:, 3] .- stressdisp[:, 3], npts, L"\sigma_{xy} \; \mathrm{(total)}")
    PyPlot.tight_layout()
    PyPlot.show()
end
ex_flamant()
