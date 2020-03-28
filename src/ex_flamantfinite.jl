using Revise
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function stylesubplots_local()
    gca().set_aspect("equal")
    PyPlot.xlim([-200, 200])
    PyPlot.ylim([-200, 200])
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
    # contourvecstress = collect(LinRange(-1e-2, 1e-2, 51))
    contourvecstress = 51

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


function ex_flamantfinite()
    obswidth = 200

    # Uniform mesh spacing
    npts = 201
    x1, y1, x2, y2 = discretizedline(-2000, 0, 2000, 0, npts)
    @show a = x1[2] - x1[1]
    mididx = Int(floor(length(x1)/2)+1)
    
    mu = 3e10
    nu = 0.25
    npts = 200
    x, y = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    nels = length(x1)
    P = 1.0
    fx = zeros(nels)
    fy = zeros(nels)
    fy[mididx] = 1.0
    
    #! Rescale forcing for smaller central element. This is strange
    fxscaled = @. fx / (2*a) #! Strange length dependent normalization
    fyscaled = @. fy / (2*a) #! Strange length dependent normalization

    #! Analytic finite Flamant solution
    # dispanalytic, stressanalytic = flamant(x, y, fx, fy, mididx, mu, nu)
    dispanalytic = zeros(length(x), 2)
    stressanalytic = zeros(length(x), 3)
    r1 = @. (x-a)^2 + y^2
    r2 = @. (x+a)^2 + y^2
    theta1 = @. atan(y, x-a)
    theta2 = @. atan(y, x+a)
    stressanalytic[:, 1] = @. -P/pi * (theta1-theta2 + y*(x-a)*r1^-2 - y*(x+a)*r2^-2)
    stressanalytic[:, 2] = @. -P/pi * (theta1-theta2 - y*(x-a)*r1^-2 + y*(x+a)*r2^-2)
    stressanalytic[:, 3] = @. -P/pi * y^2 * (r1^-2 - r2^-2)

    #! BEM solution
    els = Elements(Int(1e5))
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "line"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Direct BEM
    T, _, _ = partialsconstdispstress(slip2dispstress, els, idx["line"], idx["line"], mu, nu)
    U, _, _ = partialsconstdispstress(trac2dispstress, els, idx["line"], idx["line"], mu, nu)
    # dispall = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * interleave(fxscaled, fyscaled)
    dispall = (inv(T + 0.5 * I(size(T)[1]))) * U * interleave(fxscaled, fyscaled)

    #! Try Indirect DDM...because it worked for the Brazil test
    dispall = -U * interleave(fx, fy)

    #! Streses from tractions and induced displacements
    disptrac, stresstrac = constdispstress(trac2dispstress, x, y, els, idx["line"], fxscaled, fyscaled, mu, nu)
    dispdisp, stressdisp = constdispstress(slip2dispstress, x, y, els, idx["line"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    #! Set everything in the upper half-plane to zero because that is the exterior solution
    to_nan_idx = findall(x -> x > 0.0, y)
    dispanalytic[to_nan_idx, :] .= NaN
    stressanalytic[to_nan_idx, :] .= NaN
    disptrac[to_nan_idx, :] .= NaN
    dispdisp[to_nan_idx, :] .= NaN
    stresstrac[to_nan_idx, :] .= NaN
    stressdisp[to_nan_idx, :] .= NaN

    #! 18 panel plot
    close("all")
    # dispbem = disptrac .- dispdisp # Direct BEM
    # stressbem = stresstrac .- stressdisp # Direct BEM
    dispbem = dispdisp # DDM
    stressbem = stressdisp # DDM

    xobs = reshape(x, npts, npts)
    yobs = reshape(y, npts, npts)
    plot18_local(els, xobs, yobs, dispbem, stressbem, "BEM",
                 dispanalytic, stressanalytic, "analytic",
                 "Flamant (BEM vs. analytic)")
    tight_layout()
    show()
end
ex_flamantfinite()
