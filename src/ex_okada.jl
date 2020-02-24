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
        _, u, ∇u = ow.dc3dwrapper(0.6, [0.0, x[i], y[i]-deep], 0.0 + deep, 0.0, [-1000000, 1000000], [-0.5, 0.5], [0.0, 1.0, 0.0])
        ∇u = [∇u[2,2] ∇u[2,3] ; ∇u[3,2] ∇u[3,3]]
        ϵ = @. 0.5 * (∇u + transpose(∇u))
        σ = mu*I(2)*tr(ϵ) + 2*mu*ϵ
        dispokadass[i, 1] = u[2]
        dispokadass[i, 2] = u[3]
        stressokadass[i, 1] = σ[1, 1]
        stressokadass[i, 2] = σ[2, 2]
        stressokadass[i, 3] = σ[1, 2]

        _, u, ∇u = ow.dc3dwrapper(0.6, [0.0, x[i], y[i]-deep], 0.0 + deep, 0.0, [-1000000, 1000000], [-0.5, 0.5], [0.0, 0.0, 1.0])
        ∇u = [∇u[2,2] ∇u[2,3] ; ∇u[3,2] ∇u[3,3]]
        ϵ = @. 0.5 * (∇u + transpose(∇u))
        σ = mu*I(2)*tr(ϵ) + 2*mu*ϵ
        dispokadats[i, 1] = u[2]
        dispokadats[i, 2] = u[3]
        stressokadats[i, 1] = σ[1, 1]
        stressokadats[i, 2] = σ[2, 2]
        stressokadats[i, 3] = σ[1, 2]
    end

    # Set contour levels for displacements and stresses
    contourvecdisp = collect(LinRange(-0.5, 0.5, 51))
    contourvecstress = collect(LinRange(-1e11, 1e11, 51))
    cmap = PyPlot.get_cmap("seismic")
    
    PyPlot.close("all")
    PyPlot.figure(figsize=(30, 20))

    # BEM solutions
    PyPlot.subplot(3, 6, 1)
    PyPlot.quiver(x, y, dispbemss[:, 1], dispbemss[:, 2], units="width", color="b")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\mathbf{u} \; \mathrm{(BEM)}")

    PyPlot.subplot(3, 6, 2)
    PyPlot.contourf(xobs, yobs, reshape(dispbemss[:, 1], size(xobs)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(dispbemss[:, 1], size(xobs)), contourvecdisp, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"u_x \; \mathrm{(BEM)}")
    
    PyPlot.subplot(3, 6, 3)
    PyPlot.contourf(xobs, yobs, reshape(dispbemss[:, 2], size(xobs)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(dispbemss[:, 2], size(xobs)), contourvecdisp, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"u_y \; \mathrm{(BEM)}")

    PyPlot.subplot(3, 6, 4)
    PyPlot.contourf(xobs, yobs, reshape(stressbemss[:, 1], size(xobs)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(stressbemss[:, 1], size(xobs)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\sigma_{xx} \; \mathrm{(BEM)}")

    PyPlot.subplot(3, 6, 5)
    PyPlot.contourf(xobs, yobs, reshape(stressbemss[:, 2], size(xobs)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(stressbemss[:, 2], size(xobs)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\sigma_{yy} \; \mathrm{(BEM)}")

    PyPlot.subplot(3, 6, 6)
    PyPlot.contourf(xobs, yobs, reshape(stressbemss[:, 3], size(xobs)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(stressbemss[:, 3], size(xobs)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\sigma_{xy} \; \mathrm{(BEM)}")

    # Okada solutions
    PyPlot.subplot(3, 6, 7)
    PyPlot.quiver(x, y, dispokadass[:, 1], dispokadass[:, 2], units="width", color="b")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\mathbf{u} \; \mathrm{(Okada)}")

    PyPlot.subplot(3, 6, 8)
    PyPlot.contourf(xobs, yobs, reshape(dispbemss[:, 1], size(xobs)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(dispbemss[:, 1], size(xobs)), contourvecdisp, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"u_x \; \mathrm{(Okada)}")
    
    PyPlot.subplot(3, 6, 9)
    PyPlot.contourf(xobs, yobs, reshape(dispokadass[:, 2], size(xobs)), contourvecdisp, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(dispokadass[:, 2], size(xobs)), contourvecdisp, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"u_y \; \mathrm{(Okada)}")

    PyPlot.subplot(3, 6, 10)
    PyPlot.contourf(xobs, yobs, reshape(stressokadass[:, 1], size(xobs)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(stressokadass[:, 1], size(xobs)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\sigma_{xx} \; \mathrm{(Okada)}")

    PyPlot.subplot(3, 6, 11)
    PyPlot.contourf(xobs, yobs, reshape(stressokadass[:, 2], size(xobs)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(stressokadass[:, 2], size(xobs)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\sigma_{yy} \; \mathrm{(Okada)}")

    PyPlot.subplot(3, 6, 12)
    PyPlot.contourf(xobs, yobs, reshape(stressokadass[:, 3], size(xobs)), contourvecstress, cmap=cmap)
    PyPlot.colorbar(fraction=0.020, pad=0.05, extend="both")
    PyPlot.contour(xobs, yobs, reshape(stressokadass[:, 3], size(xobs)), contourvecstress, linewidths=0.25, colors="k")
    stylesubplots_local()
    Bem2d.plotelements(els)
    PyPlot.title(L"\sigma_{xy} \; \mathrm{(Okada)}")

    PyPlot.show()
end
ex_okada()
