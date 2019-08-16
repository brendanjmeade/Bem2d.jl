module Bem2d
using LaTeXStrings
# using Plots
using PyCall
using PyPlot

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
    a::Array{Float64, 1}
    b::Array{Float64, 1}
    sigman::Array{Float64, 1}
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
                fill(NaN, maxidx), # a::Array{Float64, 1}
                fill(NaN, maxidx), # b::Array{Float64, 1}
                fill(NaN, maxidx), # sigman::Array{Float64, 1}
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
    # display(Base.@locals())
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
        _stress[i, 1], _stress[i, 3], _, _stress[i, 2] = rotmat * stresstensor * rotmat'
    end
    return _disp, _stress
end

export plotelements
function plotelements(elements)
    for i in 1:elements.endidx
        plot([elements.x1[i], elements.x2[i]],
             [elements.y1[i], elements.y2[i]], "-k",
             linewidth=0.5)
    end
end

function stylesubplots(xlim, ylim)
    gca().set_aspect("equal")
    gca().set_xlim([xlim[1], xlim[2]])
    gca().set_ylim([ylim[1], ylim[2]])
    gca().set_xticks([xlim[1], xlim[2]])
    gca().set_yticks([ylim[1], ylim[2]])
end

function plotfields_contours(elements, xobs, yobs, idx, field, title)
    ncontours = 10
    xlim = [minimum(xobs) maximum(xobs)]
    ylim = [minimum(yobs) maximum(yobs)]
    subplot(2, 3, idx)
    scale = 5e-1
    fieldmax = maximum(@.abs(field))
    contourf(xobs, yobs, reshape(field, size(xobs)), ncontours,
        vmin=-scale * fieldmax, vmax=scale * fieldmax, cmap=plt.get_cmap("PiYG"))
    clim(-scale * fieldmax, scale * fieldmax)
    colorbar(fraction=0.020, pad=0.05, extend="both")
    contour(xobs, yobs, reshape(field, size(xobs)), ncontours,
        vmin=-scale * fieldmax, vmax=scale * fieldmax, linewidths=0.25, colors="k")
    PyPlot.title(title)
    stylesubplots(xlim, ylim)
    plotelements(elements)
end

export plotfields
function plotfields(elements, xobs, yobs, disp, stress, suptitle)
    figure(figsize = (16, 8))
    subplot(2, 3, 1)
    quiver(xobs[:], yobs[:], disp[:, 1], disp[:, 2], units="width", color="b")
    stylesubplots([minimum(xobs) maximum(xobs)], [minimum(yobs) maximum(yobs)])
    plotelements(elements)
    PyPlot.title("displacements")

    plotfields_contours(elements, xobs, yobs, 2, disp[:, 1], L"u_x")
    plotfields_contours(elements, xobs, yobs, 3, disp[:, 2], L"u_y")
    plotfields_contours(elements, xobs, yobs, 4, stress[:, 1], L"\sigma_{xx}")
    plotfields_contours(elements, xobs, yobs, 5, stress[:, 2], L"\sigma_{yy}")
    plotfields_contours(elements, xobs, yobs, 6, stress[:, 3], L"\sigma_{xy}")

    PyPlot.suptitle(suptitle)
    tight_layout()
    show()
end

export partials_constslip
function partials_constslip(elements, srcidx, obsidx, mu, nu)
    partialsdisp = zeros(2 * length(obsidx), 2 * length(srcidx))
    partialsstress = zeros(3 * length(obsidx), 2 * length(srcidx))
    partialstrac = zeros(2 * length(obsidx), 2 * length(srcidx))

    for isrc in 1:length(srcidx)
        for iobs in 1:length(obsidx)
            pd, ps, pt = zeros(2, 2), zeros(3, 2), zeros(2, 2)
            pd[:, 1], ps[:, 1] = dispstress_constslip(
                elements.xcenter[obsidx[iobs]], elements.ycenter[obsidx[iobs]],
                elements.halflength[srcidx[isrc]], mu, nu, 1, 0,
                elements.xcenter[srcidx[isrc]], elements.ycenter[srcidx[isrc]],
                elements.rotmat[srcidx[isrc], :, :], elements.rotmatinv[srcidx[isrc], :, :])
            pd[:, 2], ps[:, 2] = dispstress_constslip(
                elements.xcenter[obsidx[iobs]], elements.ycenter[obsidx[iobs]],
                elements.halflength[srcidx[isrc]], mu, nu, 0, 1,
                elements.xcenter[srcidx[isrc]], elements.ycenter[srcidx[isrc]],
                elements.rotmat[srcidx[isrc], :, :], elements.rotmatinv[srcidx[isrc], :, :])
            pt[:, 1] = stress2trac(ps[:, 1], [elements.xnormal[obsidx[iobs]] ; elements.ynormal[obsidx[iobs]]])
            pt[:, 2] = stress2trac(ps[:, 2], [elements.xnormal[obsidx[iobs]] ; elements.ynormal[obsidx[iobs]]])
            partialsdisp[2 * (iobs - 1) + 1 : 2 * (iobs - 1) + 2,
                          2 * (isrc - 1) + 1 : 2 * (isrc - 1) + 2] = pd
            partialsstress[3 * (iobs - 1) + 1 : 3 * (iobs - 1) + 3,
                            2 * (isrc - 1) + 1 : 2 * (isrc - 1) + 2] = ps
            partialstrac[2 * (iobs - 1) + 1 : 2 * (iobs - 1) + 2,
                          2 * (isrc - 1) + 1 : 2 * (isrc - 1) + 2] = pt
        end
    end
    return partialsdisp, partialsstress, partialstrac
end

end
