module Bem2d
using LaTeXStrings
using Plots

# Based on: http://julia-programming-language.2336112.n4.nabble.com/Meshgrid-function-td37003.html
export meshgrid
function meshgrid(xs, ys)
    return [xs[i] for i in 1:length(xs), j in 1:length(ys)], [ys[j] for i in 1:length(xs), j in 1:length(ys)]
end

maxidx = Int64(1e5)
println()
println("Bem2d.jl: Maximum number of elements: maxidx = ", maxidx)
println()
export Elements
mutable struct Elements
    x1::Array{Float64, 1}
    y1::Array{Float64, 1}
    x2::Array{Float64, 1}
    y2::Array{Float64, 1}
    name::Array{String, 1}
    uxconst::Array{Float64, 1}
    uyconst::Array{Float64, 1}
    uxquad::Array{Float64, 2}
    uyquad::Array{Float64, 2}
    angle::Array{Float64, 1}
    length::Array{Float64, 1}
    halflength::Array{Float64, 1}
    xcenter::Array{Float64, 1}
    ycenter::Array{Float64, 1}
    rotmat::Array{Float64, 3}
    rotmatinv::Array{Float64, 3}
    xnormal::Array{Float64, 1}
    ynormal::Array{Float64, 1}
    xnodes::Array{Float64, 2}
    ynodes::Array{Float64, 2}
    endidx::Int64
    Elements() = new(fill(NaN, maxidx), # x1
                fill(NaN, maxidx), # y1
                fill(NaN, maxidx), # x2
                fill(NaN, maxidx), # y2
                fill("", maxidx), # names
                fill(NaN, maxidx), # uxglobal::Array{Float64, 1}
                fill(NaN, maxidx), # uyglobal::Array{Float64, 1}
                fill(NaN, maxidx, 3), # uxglobalquad::Array{Float64, 2}
                fill(NaN, maxidx, 3), # uyglobalquad::Array{Float64, 2}
                fill(NaN, maxidx), # angle::Array{Float64, 1}
                fill(NaN, maxidx), # length::Array{Float64, 1}
                fill(NaN, maxidx), # halflength::Array{Float64, 1}
                fill(NaN, maxidx), # xcenter::Array{Float64, 1}
                fill(NaN, maxidx), # ycenter::Array{Float64, 1}
                fill(NaN, maxidx, 2, 2), # rotationmatrix::Array{Float64, 3}
                fill(NaN, maxidx, 2, 2), # rotationmatrixinverse::Array{Float64, 3}
                fill(NaN, maxidx), # xnormal::Array{Float64, 1}
                fill(NaN, maxidx), # ynormal::Array{Float64, 1}
                fill(NaN, maxidx, 3), # xnodes
                fill(NaN, maxidx, 3), # ynodes
                0) # endidx
end

export updateendidx!
function updateendidx!(elements)
    elements.endidx = findall(isnan, elements.x1)[1] - 1
    return nothing
end

export discretizedline
function discretizedline(xstart, ystart, xend, yend, nelements)
    npts = nelements + 1
    x = range(xstart, stop=xend, length=npts)
    y = range(ystart, stop=yend, length=npts)
    x1 = x[1:1:end-1]
    y1 = y[1:1:end-1]
    x2 = x[2:1:end]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end

# Calculate displacements and stresses for constant slip elements
export dispstress_constslip
function dispstress_constslip(x, y, a, mu, nu, xcomp, ycomp, xcenter, ycenter, rotmat, rotmatinv)
    display(Base.@locals())
    # Rotate and translate into local coordinate system with *global* slip components
    _x = zeros(length(x))
    _y = zeros(length(y))
    for i in 1:length(x)
        _x[i], _y[i] = rotmatinv * [x[i] - xcenter ; y[i] - ycenter]
    end
    _xcomp, _ycomp = rotmatinv * [xcomp ; ycomp]
    f = constantkernel(_x, _y, a, nu)
    disp, stress = slip2dispstress(_xcomp, _ycomp, f, _y, mu, nu)
    disp, stress = rotdispstress(disp, stress, rotmat)
    return disp, stress
end

# Calculate displacements and stresses for constant traction elements
export dispstress_consttrac
function dispstress_consttrac(x, y, a, mu, nu, xcomp, ycomp, xcenter, ycenter, rotmat, rotmatinv)
    # Rotate and translate into local coordinate system with *global* slip components
    _x = zeros(length(x))
    _y = zeros(length(y))
    for i in 1:length(x)
        _x[i], _y[i] = rotmat * [x[i] - xcenter ; y[i] - ycenter]
    end
    _xcomp, _ycomp = rotmatinv * [xcomp ; ycomp]
    f = constantkernel(_x, _y, a, nu)
    disp, stress = trac2dispstress(_xcomp, _ycomp, f, _y, mu, nu)
    disp, stress = rotdispstress(disp, stress, rotmatinv)
    return disp, stress
end

# Constant slip kernels from Starfield and Crouch, pages 49 and 82
function constantkernel(x, y, a, nu)
    f = zeros(length(x), 7)
    for i in 1:length(x)
        f[i, 1] = -1 / (4 * pi * (1 - nu)) * (y[i] * (atan(y[i], (x[i] - a)) - atan(y[i], (x[i] + a)))- (x[i] - a) * log(sqrt((x[i] - a)^2 + y[i]^2)) + (x[i] + a) * log(sqrt((x[i] + a)^2 + y[i]^2)))
        f[i, 2] = -1 / (4 * pi * (1 - nu)) * ((atan(y[i], (x[i] - a)) - atan(y[i], (x[i] + a))))
        f[i, 3] = 1 / (4 * pi * (1 - nu)) * (log(sqrt((x[i] - a)^2 + y[i]^2)) - log(sqrt((x[i] + a)^2 + y[i]^2)))
        f[i, 4] = 1 / (4 * pi * (1 - nu)) * (y[i] / ((x[i] - a)^2 + y[i]^2) - y[i] / ((x[i] + a)^2 + y[i]^2))
        f[i, 5] = 1 / (4 * pi * (1 - nu)) * ((x[i] - a) / ((x[i] - a)^2 + y[i]^2) - (x[i] + a) / ((x[i] + a)^2 + y[i]^2))
        f[i, 6] = 1 / (4 * pi * (1 - nu)) * (((x[i] - a)^2 - y[i]^2) / ((x[i] - a)^2 + y[i]^2)^2 - ((x[i] + a)^2 - y[i]^2) / ((x[i] + a)^2 + y[i]^2)^2)
        f[i, 7] = 2 * y[i] / (4 * pi * (1 - nu)) * ((x[i] - a) / ((x[i] - a)^2 + y[i]^2)^2 - (x[i] + a) / ((x[i] + a)^2 + y[i]^2)^2)
    end
    return f
end

# Generalization from Starfield and Crouch
export trac2dispstress
function trac2dispstress(xcomp, ycomp, f, y, mu, nu)
    disp = zeros(length(y), 2)
    stress = zeros(length(y), 3)
    _xcomp = -xcomp # For Okada consistency
    _ycomp = -ycomp # For Okada consistency
    for i in 1:length(y)
        disp[i, 1] = _xcomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f[i, 1] + y[i] * f[i, 2]) + _ycomp / (2.0 * mu) * (-y[i] * f[i, 3])
        disp[i, 2] = _xcomp / (2.0 * mu) * (-y[i] * f[i, 3]) + _ycomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f[i, 1] - y[i] * f[i, 2])
        stress[i, 1] = _xcomp * ((3.0 - 2.0 * nu) * f[i, 3] + y[i] * f[i, 4]) + _ycomp * (2.0 * nu * f[i, 2] + y[i] * f[i, 5])
        stress[i, 2] = _xcomp * (-1.0 * (1.0 - 2.0 * nu) * f[i, 3] + y[i] * f[i, 4]) + _ycomp * (2.0 * (1.0 - nu) * f[i, 2] - y[i] * f[i, 5])
        stress[i, 3] = _xcomp * (2.0 * (1.0 - nu) * f[i, 2] + y[i] * f[i, 5]) + _ycomp * ((1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4])
    end
    return disp, stress
end

# Generalization of Starfield and Crouch
export slip2dispstress
function slip2dispstress(xcomp, ycomp, f, y, mu, nu)
    disp = zeros(length(y), 2)
    stress = zeros(length(y), 3)
    _xcomp = -xcomp # For Okada consistency
    _ycomp = -ycomp # For Okada consistency

    for i in 1:length(y)
        disp[i, 1] = _xcomp * (2.0 * (1.0 - nu) * f[i, 2] - y[i] * f[i, 5]) + _ycomp * (-1.0 * (1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4])
        disp[i, 2] = _xcomp * ((1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4]) + _ycomp * (2.0 * (1 - nu) * f[i, 2] - y[i] * -f[i, 5])
        stress[i, 1] = 2.0 * _xcomp * mu * (2.0 * f[i, 4] + y[i] * f[i, 6]) + 2.0 * _ycomp * mu * (-f[i, 5] + y[i] * f[i, 7])
        stress[i, 2] = 2.0 * _xcomp * mu * (-y[i] * f[i, 6]) + 2.0 * _ycomp * mu * (-f[i, 5] - y[i] * f[i, 7])
        stress[i, 3] = 2.0 * _xcomp * mu * (-f[i, 5] + y[i] * f[i, 7]) + 2.0 * _ycomp * mu * (-y[i] * f[i, 6])
    end
    return disp, stress
end

export stress2trac
function stress2trac(stress, nvec)
    return [stress[1] stress[3] ; stress[3] stress[2]] * nvec
end

export standardize_elements!
function standardize_elements!(elements)
    updateendidx!(elements)
    for i in 1:elements.endidx
        dx = elements.x2[i] - elements.x1[i]
        dy = elements.y2[i] - elements.y1[i]
        magnitude = sqrt(dx^2 + dy^2)
        elements.angle[i] = atan(dy, dx)
        elements.length[i] = magnitude
        elements.halflength[i] = 0.5 * elements.length[i]
        elements.xcenter[i] = 0.5 * (elements.x2[i] + elements.x1[i])
        elements.ycenter[i] = 0.5 * (elements.y2[i] + elements.y1[i])
        elements.rotmat[i, :, :] = [cos(elements.angle[i]) -sin(elements.angle[i]) ; sin(elements.angle[i]) cos(elements.angle[i])]
        elements.rotmatinv[i, :, :] = [cos(-elements.angle[i]) -sin(-elements.angle[i]) ; sin(-elements.angle[i]) cos(-elements.angle[i])]
        elements.xnormal[i] = dy / magnitude
        elements.ynormal[i] = -dx / magnitude
        elements.xnodes[i, :] = [elements.xcenter[i] - (2 / 3 * dx / 2), elements.xcenter[i], elements.xcenter[i] + (2 / 3 * dx / 2)]
        elements.ynodes[i, :] = [elements.ycenter[i] - (2 / 3 * dy / 2), elements.ycenter[i], elements.ycenter[i] + (2 / 3 * dy / 2)]
    end
    return nothing
end

export rotdispstress
function rotdispstress(disp, stress, rotmat)
    # TODO: If this is slow hand expand the matrix vector multiplies
    # inplace for speed.  Some benchmarks suggest 50x speedup!
    _disp = zeros(size(disp))
    _stress = zeros(size(stress))
    for i in 1:size(stress)[1]
        _disp[i, 1], _disp[i, 2] = rotmat * [disp[i, 1] ; disp[i, 2]]
        stresstensor = [stress[i, 1] stress[i, 3] ; stress[i, 3] stress[i, 2]]
        _stress[i, 1], _stress[i, 3], _, _stress[i, 2] = rotmat' * stresstensor * rotmat
    end
    return _disp, _stress
end

function plotelements(elements)
    for i in 1:elements.endidx
        plot!([elements.x1[i]/1e3, elements.x2[i]/1e3],
              [elements.y1[i]/1e3, elements.y2[i]/1e3],
              linewidth=0.5, linecolor=:black, legend=:none)
    end
end

export plotfields
function plotfields(elements, xgrid, ygrid, disp, stress, title)
    figsize = (1000, 800)
    quiverscale = 5.0
    quiveridx = 1:length(xgrid[:])
    p1 = quiver(xgrid[quiveridx]/1e3, ygrid[quiveridx]/1e3,
                quiver=(quiverscale .* disp[quiveridx, 1], quiverscale .* disp[quiveridx, 2]),
                aspect_ratio=:equal, title=L"u" * " (" * title * ")",
                xlabel=L"x \; \mathrm{(km)}",
                ylabel=L"y \; \mathrm{(km)}", size=figsize, legend=:none,
                framestyle=:box, xtickfontsize=12, ytickfontsize=12,
                xlims=(minimum(xgrid)/1e3, maximum(xgrid)/1e3),
                ylims=(minimum(ygrid)/1e3, maximum(ygrid)/1e3))
    plotelements(elements)

    field = disp[:, 1]
    p2 = contourf(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), color=:PiYG,
                  aspect_ratio=:equal, title=L"u_x" * " (" * title * ")",
                  xlabel=L"x \; \mathrm{(km)}",
                  ylabel=L"y \; \mathrm{(km)}", size=figsize, legend=:inside,
                  framestyle=:box, xtickfontsize=12, ytickfontsize=12,
                  xlims=(minimum(xgrid)/1e3, maximum(xgrid)/1e3),
                  ylims=(minimum(ygrid)/1e3, maximum(ygrid)/1e3))
    contour!(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), linewidth = 0.5,
             linecolor=:gray)
    plotelements(elements)

    field = disp[:, 2]
    p3 = contourf(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), color=:PiYG,
                  aspect_ratio=:equal, title=L"u_y" * " (" * title * ")",
                  xlabel=L"x \; \mathrm{(km)}",
                  ylabel=L"y \; \mathrm{(km)}", size=figsize, legend=:right,
                  framestyle=:box, xtickfontsize=12, ytickfontsize=12,
                  xlims=(minimum(xgrid)/1e3, maximum(xgrid)/1e3),
                  ylims=(minimum(ygrid)/1e3, maximum(ygrid)/1e3))
    contour!(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), linewidth = 0.5,
             linecolor=:gray)
    plotelements(elements)

    field = stress[:, 1]
    p4 = contourf(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), color=:PiYG,
                  aspect_ratio=:equal, title=L"\sigma_{xx}" * " (" * title * ")",
                  xlabel=L"x \; \mathrm{(km)}",
                  ylabel=L"y \; \mathrm{(km)}", size=figsize, legend=:right,
                  framestyle=:box, xtickfontsize=12, ytickfontsize=12,
                  xlims=(minimum(xgrid)/1e3, maximum(xgrid)/1e3),
                  ylims=(minimum(ygrid)/1e3, maximum(ygrid)/1e3))
    contour!(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), linewidth = 0.5,
             linecolor=:gray)
    plotelements(elements)

    field = stress[:, 2]
    p5 = contourf(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), color=:PiYG,
                  aspect_ratio=:equal, title=L"\sigma_{yy}" * " (" * title * ")",
                  xlabel=L"x \; \mathrm{(km)}",
                  ylabel=L"y \; \mathrm{(km)}", size=figsize, legend=:right,
                  framestyle=:box, xtickfontsize=12, ytickfontsize=12,
                  xlims=(minimum(xgrid)/1e3, maximum(xgrid)/1e3),
                  ylims=(minimum(ygrid)/1e3, maximum(ygrid)/1e3))
    contour!(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), linewidth = 0.5,
             linecolor=:gray)
    plotelements(elements)

    field = stress[:, 3]
    p6 = contourf(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), color=:PiYG,
                  aspect_ratio=:equal, title=L"\sigma_{xy}" * " (" * title * ")",
                  xlabel=L"x \; \mathrm{(km)}",
                  ylabel=L"y \; \mathrm{(km)}", size=figsize, legend=:right,
                  framestyle=:box, xtickfontsize=12, ytickfontsize=12,
                  xlims=(minimum(xgrid)/1e3, maximum(xgrid)/1e3),
                  ylims=(minimum(ygrid)/1e3, maximum(ygrid)/1e3))
    contour!(xgrid/1e3, ygrid/1e3, reshape(field, size(xgrid)), linewidth = 0.5,
             linecolor=:gray)
    plotelements(elements)

    plot(p1, p2, p3, p4, p5, p6, layout=(2, 3))
    gui()
end

export partials_constslip_single
function partials_constslip_single(elements, srcidx, obsidx, mu, nu)
    pdisp = zeros(2, 2)
    pstress = zeros(3, 2)
    ptrac = zeros(2, 2)
    display(srcidx)
    display(obsidx)
    pdisp[:, 1], pstress[:, 1] = dispstress_constslip(
        elements.xcenter[obsidx], elements.ycenter[obsidx],
        elements.halflength[srcidx], mu, nu, 1, 0,
        elements.xcenter[srcidx], elements.ycenter[srcidx],
        elements.rotmat[srcidx, :, :], elements.rotmatinv[srcidx, :, :])
    pdisp[:, 2], pstress[:, 2] = dispstress_constslip(
        elements.xcenter[obsidx], elements.ycenter[obsidx],
        elements.halflength[srcidx], mu, nu, 0, 1,
        elements.xcenter[srcidx], elements.ycenter[srcidx],
        elements.rotmat[srcidx, :, :], elements.rotmatinv[srcidx, :, :])
    display(pdisp)

    ptrac[:, 1] = stress2trac(pstress[:, 1], [elements.xnormal[obsidx] ; elements.ynormal[obsidx]])
    ptrac[:, 2] = stress2trac(pstress[:, 2], [elements.xnormal[obsidx] ; elements.ynormal[obsidx]])
    return pdisp, pstress, ptrac
end

export partials_constslip
function partials_constslip(elements, srcidx, obsidx, mu, nu)
    partials_disp = zeros(2 * length(obsidx), 2 * length(srcidx))
    partials_stress = zeros(3 * length(obsidx), 2 * length(srcidx))
    partials_trac = zeros(2 * length(obsidx), 2 * length(srcidx))

    for isrc in 1:length(srcidx)
        for iobs in 1:length(obsidx)
            # TODO: Should I just move partials_constslip_single into here???
            pd, ps, pt = partials_constslip_single(elements, srcidx[isrc], obsidx[iobs], mu, nu)
            partials_disp[2 * (iobs - 1) + 1 : 2 * (iobs - 1) + 2,
                          2 * (isrc - 1) + 1 : 2 * (isrc - 1) + 2] = pd
            partials_stress[3 * (iobs - 1) + 1 : 3 * (iobs - 1) + 3,
                            2 * (isrc - 1) + 1 : 2 * (isrc - 1) + 2] = ps
            partials_trac[2 * (iobs - 1) + 1 : 2 * (iobs - 1) + 2,
                          2 * (isrc - 1) + 1 : 2 * (isrc - 1) + 2] = pt
        end
    end
    return partials_disp, partials_stress, partials_trac
end

end
