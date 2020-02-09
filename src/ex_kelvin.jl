using Revise
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function local_subplot(x, y, mat, npts, title_string)
    fontsize = 20
    contour_levels = 10
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

    # Observation coordinates
    npts = 30
    obswidth = 1000
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)
    fx = 0.0
    fy = 1.0

    # Try the Kelvin solution from Crouch and Starfield section 4.2 (z-line load in full space)
    C = 1/(4*pi*(1-nu))
    r = @. sqrt(x^2+y^2)
    g = @. -C * log(r)
    gx = @. -C * x/(x^2+y^2)
    gy = @. -C * y/(x^2+y^2)
    gxy = @. C * 2*x*y/(x^2+y^2)^2
    gxx = @. C * (x^2-y^2)/(x^2+y^2)^2
    gyy = -gxx
    ux_kelvin = @. fx/(2*mu)*((3-4*nu)*g-x*gx) + fy/(2*mu)*(-y*gx)
    uy_kelvin = @. fx/(2*mu)*(-x*gy) + fy/(2*mu)*((3-4*nu)*g-y*gy)
    σxx_kelvin = @. fy * (2*nu*gy-y*gxx) # Need to add fx terms
    σyy_kelvin = @. fy * (2*(1-nu)*gy-y*gyy)
    σxy_kelvin = @. fy * ((1-2*nu)*gx-y*gxy)

    # BEM solution
    els = Bem2d.Elements(Int(2))
    els.x1[1] = -1e-10
    els.y1[1] = 0.0
    els.x2[1] = 1e-10
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
    disptrac, stresstrac = Bem2d.constdispstress(trac2dispstress, x, y, els, idx["point"], fx, fy, mu, nu)

    # Stresses from traction induced displacements
    dispdisp, stressdisp = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["point"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    fontsize = 20
    PyPlot.close("all")
    PyPlot.figure(figsize=(40,20))
    # Analytic Kelvin
    PyPlot.subplot(3, 6, 1)
    quiver(x[:], y[:], ux_kelvin, uy_kelvin, units = "width", color = "b")
    PyPlot.title("displacement vectors", fontsize=fontsize)
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])
    PyPlot.xticks([])
    PyPlot.yticks([])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)

    PyPlot.subplot(3, 6, 2)
    local_subplot(x, y, ux_kelvin, npts, L"u_x \; (Kelvin, CS)")
    PyPlot.subplot(3, 6, 3)
    local_subplot(x, y, uy_kelvin, npts, L"u_y \; (Kelvin, CS)")
    PyPlot.subplot(3, 6, 7)
    local_subplot(x, y, σxx_kelvin, npts, L"\sigma_{xx} (Kelvin, CS)")
    PyPlot.subplot(3, 6, 8)
    local_subplot(x, y, σyy_kelvin, npts, L"\sigma_{yy} (Kelvin, CS)")
    PyPlot.subplot(3, 6, 9)
    local_subplot(x, y, σxy_kelvin, npts, L"\sigma_{xy} (Kelvin, CS)")

    # BEM Kelvin
    PyPlot.subplot(3, 6, 4)
    quiver(x[:], y[:], disptrac[:, 1], disptrac[:, 2], units = "width", color = "b")
    PyPlot.title("displacement vectors", fontsize=fontsize)
    PyPlot.xlim([-1000, 1000])
    PyPlot.ylim([-1000, 1000])
    PyPlot.xticks([])
    PyPlot.yticks([])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
    PyPlot.subplot(3, 6, 5)
    local_subplot(x, y, disptrac[:, 1], npts, L"u_x \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 6)
    local_subplot(x, y, disptrac[:, 2], npts, L"u_y \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 10)
    local_subplot(x, y, stresstrac[:, 1], npts, L"\sigma_{xx} \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 11)
    local_subplot(x, y, stresstrac[:, 2], npts, L"\sigma_{yy} \; \mathrm{(traction \; BEM)}")
    PyPlot.subplot(3, 6, 12)
    local_subplot(x, y, stressdisp[:, 3], npts, L"\sigma_{xy} \; \mathrm{(traction \; BEM)}")
    # PyPlot.subplot(3, 6, 16)
    # local_subplot(x, y, stresstrac[:, 1]+stressdisp[:, 1], npts, L"\sigma_{xx} \; \mathrm{(total)}")
    # PyPlot.subplot(3, 6, 17)
    # local_subplot(x, y, stresstrac[:, 2]+stressdisp[:, 2], npts, L"\sigma_{yy} \; \mathrm{(total)}")
    # PyPlot.subplot(3, 6, 18)
    # local_subplot(x, y, stresstrac[:, 3]+stressdisp[:, 3], npts, L"\sigma_{xy} \; \mathrm{(total)}")

    PyPlot.tight_layout()
    PyPlot.show()

end
ex_flamant()
