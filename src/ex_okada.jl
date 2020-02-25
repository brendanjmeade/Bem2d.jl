using Revise
using PyCall
using PyPlot
using LinearAlgebra
using Bem2d

function stylesubplots_local()
    gca().set_aspect("equal")
    gca().set_xticks([])
    gca().set_yticks([])
    return nothing
end

function plot18(els, x, y, disp1, stress1, string1, disp2, stress2, string2, title_string)
    # Set contour levels for displacements and stresses
    contourvecdisp = collect(LinRange(-0.5, 0.5, 51))
    contourvecstress = collect(LinRange(-1e11, 1e11, 51))
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
    PyPlot.contourf(x, y, reshape(disp1[:, 1], size(x)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 1], size(x)), contourvecdisp, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("ux (" * string1 * ")", fontsize=fontsize)
    
    PyPlot.subplot(3, 6, 3)
    PyPlot.contourf(x, y, reshape(disp1[:, 2], size(x)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 2], size(x)), contourvecdisp, linewidths=0.25, colors="k")
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
    PyPlot.contourf(x, y, reshape(disp2[:, 1], size(x)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp2[:, 1], size(x)), contourvecdisp, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("ux (" * string2 * ")", fontsize=fontsize)
    
    PyPlot.subplot(3, 6, 9)
    PyPlot.contourf(x, y, reshape(disp2[:, 2], size(x)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp2[:, 2], size(x)), contourvecdisp, linewidths=0.25, colors="k")
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
    PyPlot.contourf(x, y, reshape(disp1[:, 1]-disp2[:, 1], size(x)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 1]-disp2[:, 1], size(x)), contourvecdisp, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title("uy (" * string1 * "-" * string2 * ")", fontsize=fontsize)
    
    PyPlot.subplot(3, 6, 15)
    PyPlot.contourf(x, y, reshape(disp1[:, 2]-disp2[:, 2], size(x)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(x, y, reshape(disp1[:, 2]-disp2[:, 2], size(x)), contourvecdisp, linewidths=0.25, colors="k")
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


function ex_okada()
    mu = 30e9
    nu = 0.25

    # Flat fault
    nfault = 1
    x1, y1, x2, y2 = Bem2d.discretizedline(-0.5, 0, 0.5, 0, nfault)
    els = Bem2d.Elements(Int(1e5))
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    obswidth = 3
    npts = 99
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    xobs = reshape(x, npts, npts)
    yobs = reshape(y, npts, npts)
    
    # Constant slip fault
    dispbemss, stressbemss = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["fault"], 1, 0, mu, nu)
    dispbemts, stressbemts = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["fault"], 0, 1, mu, nu)

    # Okada solution
    ow = pyimport("okada_wrapper") # from okada_wrapper import dc3dwrapper
    dispokadass = zeros(length(x), 2)
    stressokadass = zeros(length(x), 3)
    dispokadats = zeros(length(x), 2)
    stressokadats = zeros(length(x), 3)
    deep = 10000.0
    for i in 1:length(x)
        # Strike-slip
        _, u, gradu = ow.dc3dwrapper((2.0/3.0), [0.0, x[i], y[i]-deep], 0.0 + deep, 0.0, [-1000000, 1000000], [-0.5, 0.5], [0.0, 1.0, 0.0])
        strain = @. 0.5 * (gradu' + gradu)
        stress = mu*LinearAlgebra.I(3)*tr(strain) + 2.0*mu*strain
        dispokadass[i, 1] = u[2]
        dispokadass[i, 2] = u[3]
        stressokadass[i, 1] = stress[2, 2]
        stressokadass[i, 2] = stress[3, 3]
        stressokadass[i, 3] = stress[2, 3]

        # Tensile-slip
        _, u, gradu = ow.dc3dwrapper((2.0/3.0), [0.0, x[i], y[i]-deep], 0.0 + deep, 0.0, [-1000000, 1000000], [-0.5, 0.5], [0.0, 0.0, 1.0])
        strain = @. 0.5 * (gradu' + gradu)
        stress = mu*LinearAlgebra.I(3)*tr(strain) + 2.0*mu*strain
        dispokadats[i, 1] = u[2]
        dispokadats[i, 2] = u[3]
        stressokadats[i, 1] = stress[2, 2]
        stressokadats[i, 2] = stress[3, 3]
        stressokadats[i, 3] = stress[2, 3]
    end
    
    PyPlot.close("all")
    plot18(els, xobs, yobs, dispbemss, stressbemss, "BEM", dispokadass, stressokadass, "Okada", "strike-slip (BEM vs. Okada)")
    plot18(els, xobs, yobs, dispbemts, stressbemts, "BEM", dispokadats, stressokadats, "Okada", "tensile-slip (BEM vs. Okada)")
end
ex_okada()
