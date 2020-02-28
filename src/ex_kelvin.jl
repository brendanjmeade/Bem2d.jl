using Revise
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d


function stylesubplots_local()
    gca().set_aspect("equal")
    gca().set_xticks([])
    gca().set_yticks([])
    return nothing
end

function plot18(els, x, y, disp1, stress1, string1, disp2, stress2, string2, title_string)
    # Set contour levels for displacements and stresses
    contourvecdisp = collect(LinRange(-5e-10, 5e-10, 51))
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

function ex_kelvin()
    mu = 3e10
    nu = 0.25

    # Observation coordinates
    npts = 200
    obswidth = 1000
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    tracx = 1.0
    tracy = 0.0

    # Try the Kelvin solution from Crouch and Starfield section 4.2 (z-line load in full space)
    C = 1/(4*pi*(1-nu))
    r = @. sqrt(x^2+y^2)
    g = @. -C * log(r)
    gx = @. -C * x/(x^2+y^2)
    gy = @. -C * y/(x^2+y^2)
    gxy = @. C * 2*x*y/(x^2+y^2)^2
    gxx = @. C * (x^2-y^2)/(x^2+y^2)^2
    gyy = -gxx
    ux_kelvin = @. tracx/(2*mu)*((3-4*nu)*g-x*gx) + tracy/(2*mu)*(-y*gx)
    uy_kelvin = @. tracx/(2*mu)*(-x*gy) + tracy/(2*mu)*((3-4*nu)*g-y*gy)
    σxx_kelvin = @. tracx*(2*(1-nu)*gx-x*gxx) + tracy*(2*nu*gy-y*gxx)
    σyy_kelvin = @. tracx*(2*nu*gx-x*gyy) + tracy*(2*(1-nu)*gy-y*gyy)
    σxy_kelvin = @. tracx*((1-2*nu)*gy-x*gxy) + tracy*((1-2*nu)*gx-y*gxy)


    # Try the finite length "Kelvin" solution from Crouch and Starfied, section 4.3
    a = 0.5
    f = @. -C * (y*(atan(y, x-a)-atan(y, x+a)) - (x-a)*log(sqrt((x-a)^2 +y^2)) + (x+a)*log(sqrt((x+a)^2 +y^2)))
    fx = @. C * (log(sqrt((x-a)^2 +y^2)) - log(sqrt((x+a)^2 +y^2)))
    fy = @. -C * (atan(y,x-a) - atan(y,x+a))
    fxy = @. C * (y/((x-a)^2+y^2) - y/((x+a)^2+y^2))
    fxx = @. C * ((x-a)/((x-a)^2+y^2) - (x+a)/((x+a)^2+y^2))
    fyy = -fxx
    ux_segment = @. tracx/(2*mu)*((3-4*nu)*f-y*fx) + tracy/(2*mu)*(-y*fx)
    uy_segment = @. tracx/(2*mu)*(-y*fy) + tracy/(2*mu)*((3-4*nu)*f-y*fy)
    σxx_segment = @. tracx*((3-2*nu)*fx+y*fxy) + tracy*(2*nu*fy+y*fyy)
    σyy_segment = @. tracx*(-(1-2*nu)*fx-y*fxy) + tracy*(2*(1-nu)*fy-f*fyy)
    σxy_segment = @. tracx*(2*(1-nu)*fy+y*fyy ) + tracy*((1-2*nu)fx-y*fxy)

    # BEM solution
    els = Bem2d.Elements(Int(2))
    els.x1[1] = -a
    els.y1[1] = 0.0
    els.x2[1] = a
    els.y2[1] = 0.0
    els.name[1] = "point"
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    # Streses from tractions
    disptrac, stresstrac = Bem2d.constdispstress(trac2dispstress, x, y, els, idx["point"], tracx, tracy, mu, nu)

    fontsize = 20
    PyPlot.close("all")

    # 18 panel plot
    dispkelvin = [ux_kelvin[:] uy_kelvin[:]]
    stresskelvin = [σxx_kelvin[:] σyy_kelvin[:] σxy_kelvin[:]]
    xobs = reshape(x, npts, npts)
    yobs = reshape(y, npts, npts)
    plot18(els, xobs, yobs, disptrac, stresstrac, "BEM", dispkelvin, stresskelvin, "Kelvin point", "fx (BEM vs. Kelvin point)")

    PyPlot.figure(figsize=(40,20))
    # Analytic Kelvin
    PyPlot.subplot(3, 6, 1)
    quiver(x[:], y[:], ux_kelvin, uy_kelvin, units = "width", color = "b")
    PyPlot.title(L"\mathbf{u} \; \mathrm{(CS \; Kelvin)}", fontsize=fontsize)
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])
    PyPlot.xticks([])
    PyPlot.yticks([])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
    PyPlot.subplot(3, 6, 2)
    local_subplot(x, y, ux_kelvin, npts, L"u_x \; \mathrm{(CS \; Kelvin)}")
    PyPlot.subplot(3, 6, 3)
    local_subplot(x, y, uy_kelvin, npts, L"u_y \; \mathrm{(CS \; Kelvin)}")
    PyPlot.subplot(3, 6, 4)
    local_subplot(x, y, σxx_kelvin, npts, L"\sigma_{xx} \; \mathrm{(CS \; Kelvin)}")
    PyPlot.subplot(3, 6, 5)
    local_subplot(x, y, σyy_kelvin, npts, L"\sigma_{yy} \; \mathrm{(CS \; Kelvin)}")
    PyPlot.subplot(3, 6, 6)
    local_subplot(x, y, σxy_kelvin, npts, L"\sigma_{xy} \; \mathrm{(CS \; Kelvin)}")

    # BEM Kelvin
    PyPlot.subplot(3, 6, 7)
    quiver(x[:], y[:], disptrac[:, 1], disptrac[:, 2], units = "width", color = "b")
    PyPlot.title(L"\mathbf{u} \; \mathrm{(traction \; BEM)}", fontsize=fontsize)
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])
    PyPlot.xticks([])
    PyPlot.yticks([])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
    PyPlot.subplot(3, 6, 8)
    local_subplot(x, y, disptrac[:, 1], npts, L"u_x \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 9)
    local_subplot(x, y, disptrac[:, 2], npts, L"u_y \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 10)
    local_subplot(x, y, stresstrac[:, 1], npts, L"\sigma_{xx} \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 11)
    local_subplot(x, y, stresstrac[:, 2], npts, L"\sigma_{yy} \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 12)
    local_subplot(x, y, stresstrac[:, 3], npts, L"\sigma_{xy} \; \mathrm{(traction \; BEM)}")

    # Residuals
    PyPlot.subplot(3, 6, 13)
    quiver(x[:], y[:], ux_kelvin-disptrac[:,1], uy_kelvin-disptrac[:,2], units = "width", color = "b")
    PyPlot.title(L"\mathbf{u} \; \mathrm{(residual)}", fontsize=fontsize)
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])
    PyPlot.xticks([])
    PyPlot.yticks([])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
    PyPlot.subplot(3, 6, 14)
    local_subplot(x, y, ux_kelvin-disptrac[:,1], npts, L"u_x \; (\mathrm{residual})")
    PyPlot.subplot(3, 6, 15)
    local_subplot(x, y, uy_kelvin-disptrac[:,2], npts, L"u_y \; (\mathrm{residual})")
    PyPlot.subplot(3, 6, 16)
    local_subplot(x, y, σxx_kelvin-stresstrac[:,1], npts, L"\sigma_{xx} \; (\mathrm{residual})")
    PyPlot.subplot(3, 6, 17)
    local_subplot(x, y, σyy_kelvin-stresstrac[:,2], npts, L"\sigma_{yy} \; (\mathrm{residual})")
    PyPlot.subplot(3, 6, 18)
    local_subplot(x, y, σxy_kelvin-stresstrac[:,3], npts, L"\sigma_{xy} \; (\mathrm{residual})")

    # "Kelvin" line segment solution
    # PyPlot.subplot(3, 6, 13)
    # quiver(x[:], y[:], ux_segment, uy_segment, units = "width", color = "b")
    # PyPlot.title(L"\mathbf{u} \; \mathrm{(CS \; segment)}", fontsize=fontsize)
    # PyPlot.xlim([-1000, 1000])
    # PyPlot.ylim([-1000, 1000])
    # PyPlot.xticks([])
    # PyPlot.yticks([])
    # PyPlot.gca().set_aspect("equal")
    # PyPlot.gca().tick_params(labelsize=fontsize)
    # PyPlot.subplot(3, 6, 14)
    # local_subplot(x, y, ux_segment, npts, L"u_x \; \mathrm{(CS \; segment)}")
    # PyPlot.subplot(3, 6, 15)
    # local_subplot(x, y, uy_segment, npts, L"u_y \; \mathrm{(CS \; segment)}")
    # PyPlot.subplot(3, 6, 16)
    # local_subplot(x, y, σxx_segment, npts, L"\sigma_{xx} \; \mathrm{(CS \; segment)}")
    # PyPlot.subplot(3, 6, 17)
    # local_subplot(x, y, σyy_segment, npts, L"\sigma_{xx} \; \mathrm{(CS \; segment)}")
    # PyPlot.subplot(3, 6, 18)
    # local_subplot(x, y, σxy_segment, npts, L"\sigma_{xx} \; \mathrm{(CS \; segment)}")

    PyPlot.tight_layout()
    PyPlot.show()

end
ex_kelvin()
