using Revise
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function stylesubplots_local()
    gca().set_aspect("equal")
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])
    # gca().# set_xticks([])
    # gca().
    # set_yticks([])
    return nothing
end

function plot18_local(els, x, y, disp1, stress1, string1, disp2, stress2, string2, title_string)    
    # Set contour levels for displacements and stresses
    contourvecdispx = collect(LinRange(-1e-11, 1e-11, 51))
    # contourvecdispy = collect(LinRange(-1e-11, 1e-11, 51))
    contourvecdispy = collect(LinRange(-1e-10, 1e-10, 51))
    # contourvecdispy = 51
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
    PyPlot.tight_layout()
    PyPlot.show()
end

#! Analytic Flamant solution
function flamant(x, y, fx, fy, mididx, mu, nu)
    r = @. sqrt(x^2 + y^2)
    θ = @. rad2deg(atan(y, x))

    #! Displacement solutions (https://en.wikipedia.org/wiki/Flamant_solution)
    κplanestrain = 3-4*nu
    κplanestress = (3-nu)/(1+nu)
    κ = κplanestrain
    ur = @. -1/(4*pi*mu) * 
         (fx[mididx]*((κ-1)*deg2rad(θ)*sind(θ) - cosd(θ) + (κ+1)*cosd(θ)*log(r)) +
          fy[mididx]*((κ-1)*deg2rad(θ)*cosd(θ) + sind(θ) - (κ+1)*sind(θ)*log(r)))
    uθ = @. -1/(4*pi*mu) * 
         (fx[mididx]*((κ-1)*deg2rad(θ)*cosd(θ) - sind(θ) - (κ+1)*sind(θ)*log(r)) +
          fy[mididx]*((κ-1)*deg2rad(θ)*sind(θ) + cosd(θ) + (κ+1)*cosd(θ)*log(r)))
 

    #! Convert from cylindrical to Cartesian coordinates
    dispanalytic = zeros(length(x), 2)
    for i in 1:length(x)
        dispanalytic[i, 1] = ur[i] * cos(uθ[i])
        dispanalytic[i, 2] = ur[i] * sin(uθ[i])
    end
    dispanalytic = zeros(length(x), 2)

    #! Analytic solution from Wikipedia in cylindrical coordinates
    σrr = @. -2.0/(pi*r) * (fx[mididx]*cosd(θ) + fy[mididx]*sind(θ))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))

    #! Convert to from cylindrical to Cartesian coordinates
    σxx_flamant = zeros(length(x))
    σyy_flamant = zeros(length(x))
    σxy_flamant = zeros(length(x))
    for i in 1:length(x) # Project a single matrix to Cartesian coordinates
        cylindrical_stress_tensor = [σrr[i] σrθ[i] ; σrθ[i] σθθ[i]]
        transformation_matrix = [cosd(θ[i]) -sind(θ[i]) ; sind(θ[i]) cosd(θ[i])]
        cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transformation_matrix'
        σxx_flamant[i] = cartesian_stress_tensor[1, 1]
        σyy_flamant[i] = cartesian_stress_tensor[2, 2]
        σxy_flamant[i] = cartesian_stress_tensor[1, 2]
    end

    stressanalytic = [σxx_flamant σyy_flamant σxy_flamant]
    return dispanalytic, stressanalytic
end

function ex_flamant()
    obswidth = 1000

    #! Try wierd symmetric quadratic spacing with small elements i the middle
    a = 2.00
    nquad = 2000
    x1 = zeros(3 + 2 * (nquad-1))
    x2 = zeros(3 + 2 * (nquad-1))
    quadspacing = collect(LinRange(1, sqrt(30*obswidth), nquad).^2)
    x1[1:1:nquad-1] = -reverse(quadspacing[2:1:end], dims=1)
    x1[nquad:1:nquad+2] = [-1.0 -a a]
    x1[nquad+3:1:end] = quadspacing[1:1:end-1]
    x2[1:1:nquad-1] = -reverse(quadspacing[1:1:end-1], dims=1)
    x2[nquad:1:nquad+2] = [-a a 1.0]
    x2[nquad+3:1:end] = quadspacing[2:1:end]
    y1 = zeros(size(x1))
    y2 = zeros(size(x1))
    mididx = Int(floor(length(x1)/2)+1)

    # Uniform mesh spacing
    # npts = 2001
    # x1, y1, x2, y2 = discretizedline(-2000, 0, 2000, 0, npts)
    # mididx = Int(floor(length(x1)/2)+1)
    
    mu = 3e10
    nu = 0.25
    npts = 50
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    nels = length(x1)
    fx = zeros(nels)
    fy = zeros(nels)
    fy[mididx] = 1.0
    
    # Rescale forcing for smaller central element.
    # This is strange
    fxscaled = fx / (2*a) #! Strange length dependent normalization
    fyscaled = fy / (2*a) #! Strange length dependent normalization

    #! Analytic Flamant solution
    dispanalytic, stressanalytic = flamant(x, y, fx, fy, mididx, mu, nu)

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
    # xdisp = zeros(els.endidx)
    # ydisp = zeros(els.endidx)
    # for i in 1:els.endidx # Calcuate the x and y components of the tractions
    #     T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["line"][i], idx["line"][i], mu, nu)
    #     U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["line"][i], idx["line"][i], mu, nu)
    #     xdisp[i], ydisp[i] = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * [fx[i]; fy[i]]
    # end
    # displocal = Bem2d.interleave(xdisp, ydisp)

    # Given the other tractions calcuate the induced displacements on the boudaries
    T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["line"], idx["line"], mu, nu)
    U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["line"], idx["line"], mu, nu)
    dispall = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * Bem2d.interleave(fxscaled, fyscaled)

    # Streses from tractions and induced displacements
    disptrac, stresstrac = Bem2d.constdispstress(trac2dispstress, x, y, els, idx["line"], fxscaled, fyscaled, mu, nu)
    dispdisp, stressdisp = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["line"], dispall[1:2:end], dispall[2:2:end], mu, nu)
    
    #! Set everything in the upper half-plane to zero because that is the exterior solution
    to_nan_idx = findall(x -> x > 0.0, y)
    dispanalytic[to_nan_idx, :] .= NaN
    stressanalytic[to_nan_idx, :] .= NaN
    disptrac[to_nan_idx, :] .= NaN
    dispdisp[to_nan_idx, :] .= NaN
    stresstrac[to_nan_idx, :] .= NaN
    stressdisp[to_nan_idx, :] .= NaN

    #! 18 panel plot
    PyPlot.close("all")
    dispbem = disptrac .- dispdisp
    stressbem = stresstrac .- stressdisp
    xobs = reshape(x, npts, npts)
    yobs = reshape(y, npts, npts)
    plot18_local(els, xobs, yobs, dispbem, stressbem, "BEM", dispanalytic, stressanalytic, "analytic", "Flamant (BEM vs. analytic)")
    PyPlot.tight_layout()
    PyPlot.show()
end
ex_flamant()
