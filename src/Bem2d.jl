module Bem2d
using LaTeXStrings
using PyCall
using PyPlot
include("Bem2dQuadKernels.jl")
using Bem2dQuadKernels

# Based on: http://julia-programming-language.2336112.n4.nabble.com/Meshgrid-function-td37003.html
export meshgrid
function meshgrid(xs, ys)
    return [xs[i] for i in 1:length(xs), j in 1:length(ys)], [ys[j] for i in 1:length(xs), j in 1:length(ys)]
end

maxidx = Int64(1e5)
println()
println("Bem2d.jl: Maxiμm number of elements: maxidx = ", maxidx)
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
export constslip
function constslip(x, y, a, μ, ν, xcomp, ycomp, xcenter, ycenter, rotmat, rotmatinv)
    # Rotate and translate into local coordinate system with *global* slip components
    _x = zeros(length(x))
    _y = zeros(length(y))
    for i in 1:length(x)
        _x[i], _y[i] = rotmatinv * [x[i] - xcenter ; y[i] - ycenter]
    end
    _xcomp, _ycomp = rotmatinv * [xcomp ; ycomp]
    f = constkernel(_x, _y, a, ν)
    u, σ = slip2uσ(_xcomp, _ycomp, f, _y, μ, ν)
    u, σ = rotuσ(u, σ, rotmat)
    return u, σ
end

# Calculate displacements and stresses for constant traction elements
export consttrac
function consttrac(x, y, a, μ, ν, xcomp, ycomp, xcenter, ycenter, rotmat, rotmatinv)
    # Rotate and translate into local coordinate system with *global* slip components
    _x = zeros(length(x))
    _y = zeros(length(y))
    for i in 1:length(x)
        _x[i], _y[i] = rotmat * [x[i] - xcenter ; y[i] - ycenter]
    end
    _xcomp, _ycomp = rotmatinv * [xcomp ; ycomp]
    f = constkernel(_x, _y, a, ν)
    u, σ = trac2uσ(_xcomp, _ycomp, f, _y, μ, ν)
    u, σ = rotuσ(u, σ, rotmatinv)
    return u, σ
end

# Calculate displacements and stresses for constant slip elements
export quadslip
function quadslip(x, y, a, μ, ν, xcomp, ycomp, xcenter, ycenter, rotmat, rotmatinv)
    # Rotate and translate into local coordinate system with *global* slip components
    _x = zeros(length(x))
    _y = zeros(length(y))
    for i in 1:length(x)
        _x[i], _y[i] = rotmatinv * [x[i] - xcenter ; y[i] - ycenter]
    end
    f = quadkernel_farfield(x, y, a, ν)
    u = zeros(length(x), 2)
    σ = zeros(length(x), 2)
    for i in 1:3
        # f = f_all[:, i, :] # TODO: How to update indexing from python???
        _xcomp, _ycomp = rotmatinv * [xcomp[i] ; ycomp[i]]
        _u, _σ = slip2uσ(_xcomp, _ycomp, f[:, i, :], _y, μ, ν)
        _u, _σ = rotuσ(_u, _σ, rotmat)
        u += _u
        σ += _σ
    end
    return u, σ
end


# Constant slip kernels from Starfield and Crouch, pages 49 and 82
function constkernel(x, y, a, ν)
    f = zeros(length(x), 7)
    for i in 1:length(x)
        f[i, 1] = -1 / (4 * π * (1 - ν)) * (y[i] * (atan(y[i], (x[i] - a)) - atan(y[i], (x[i] + a)))- (x[i] - a) * log(sqrt((x[i] - a)^2 + y[i]^2)) + (x[i] + a) * log(sqrt((x[i] + a)^2 + y[i]^2)))
        f[i, 2] = -1 / (4 * π * (1 - ν)) * ((atan(y[i], (x[i] - a)) - atan(y[i], (x[i] + a))))
        f[i, 3] = 1 / (4 * π * (1 - ν)) * (log(sqrt((x[i] - a)^2 + y[i]^2)) - log(sqrt((x[i] + a)^2 + y[i]^2)))
        f[i, 4] = 1 / (4 * π * (1 - ν)) * (y[i] / ((x[i] - a)^2 + y[i]^2) - y[i] / ((x[i] + a)^2 + y[i]^2))
        f[i, 5] = 1 / (4 * π * (1 - ν)) * ((x[i] - a) / ((x[i] - a)^2 + y[i]^2) - (x[i] + a) / ((x[i] + a)^2 + y[i]^2))
        f[i, 6] = 1 / (4 * π * (1 - ν)) * (((x[i] - a)^2 - y[i]^2) / ((x[i] - a)^2 + y[i]^2)^2 - ((x[i] + a)^2 - y[i]^2) / ((x[i] + a)^2 + y[i]^2)^2)
        f[i, 7] = 2 * y[i] / (4 * π * (1 - ν)) * ((x[i] - a) / ((x[i] - a)^2 + y[i]^2)^2 - (x[i] + a) / ((x[i] + a)^2 + y[i]^2)^2)
    end
    return f
end


# Generalization from Starfield and Crouch
export trac2uσ
function trac2uσ(xcomp, ycomp, f, y, μ, ν)
    u = zeros(length(y), 2)
    σ = zeros(length(y), 3)
    _xcomp = -xcomp # For Okada consistency
    _ycomp = -ycomp # For Okada consistency
    for i in 1:length(y)
        u[i, 1] = _xcomp / (2.0 * μ) * ((3.0 - 4.0 * ν) * f[i, 1] + y[i] * f[i, 2]) + _ycomp / (2.0 * μ) * (-y[i] * f[i, 3])
        u[i, 2] = _xcomp / (2.0 * μ) * (-y[i] * f[i, 3]) + _ycomp / (2.0 * μ) * ((3.0 - 4.0 * ν) * f[i, 1] - y[i] * f[i, 2])
        σ[i, 1] = _xcomp * ((3.0 - 2.0 * ν) * f[i, 3] + y[i] * f[i, 4]) + _ycomp * (2.0 * ν * f[i, 2] + y[i] * f[i, 5])
        σ[i, 2] = _xcomp * (-1.0 * (1.0 - 2.0 * ν) * f[i, 3] + y[i] * f[i, 4]) + _ycomp * (2.0 * (1.0 - ν) * f[i, 2] - y[i] * f[i, 5])
        σ[i, 3] = _xcomp * (2.0 * (1.0 - ν) * f[i, 2] + y[i] * f[i, 5]) + _ycomp * ((1.0 - 2.0 * ν) * f[i, 3] - y[i] * f[i, 4])
    end
    return u, σ
end

# Generalization of Starfield and Crouch
export slip2uσ
function slip2uσ(xcomp, ycomp, f, y, μ, ν)
    u = zeros(length(y), 2)
    σ = zeros(length(y), 3)
    _xcomp = -xcomp # For Okada consistency
    _ycomp = -ycomp # For Okada consistency

    for i in 1:length(y)
        u[i, 1] = _xcomp * (2.0 * (1.0 - ν) * f[i, 2] - y[i] * f[i, 5]) + _ycomp * (-1.0 * (1.0 - 2.0 * ν) * f[i, 3] - y[i] * f[i, 4])
        u[i, 2] = _xcomp * ((1.0 - 2.0 * ν) * f[i, 3] - y[i] * f[i, 4]) + _ycomp * (2.0 * (1 - ν) * f[i, 2] - y[i] * -f[i, 5])
        σ[i, 1] = 2.0 * _xcomp * μ * (2.0 * f[i, 4] + y[i] * f[i, 6]) + 2.0 * _ycomp * μ * (-f[i, 5] + y[i] * f[i, 7])
        σ[i, 2] = 2.0 * _xcomp * μ * (-y[i] * f[i, 6]) + 2.0 * _ycomp * μ * (-f[i, 5] - y[i] * f[i, 7])
        σ[i, 3] = 2.0 * _xcomp * μ * (-f[i, 5] + y[i] * f[i, 7]) + 2.0 * _ycomp * μ * (-y[i] * f[i, 6])
    end
    return u, σ
end


export σ2t
function σ2t(stress, nvec)
    return [stress[1] stress[3] ; stress[3] stress[2]] * nvec
end

export standardize_elements!
function standardize_elements!(els)
    updateendidx!(els)
    for i in 1:els.endidx
        dx = els.x2[i] - els.x1[i]
        dy = els.y2[i] - els.y1[i]
        magnitude = sqrt(dx^2 + dy^2)
        els.angle[i] = atan(dy, dx)
        els.length[i] = magnitude
        els.halflength[i] = 0.5 * els.length[i]
        els.xcenter[i] = 0.5 * (els.x2[i] + els.x1[i])
        els.ycenter[i] = 0.5 * (els.y2[i] + els.y1[i])
        els.rotmat[i, :, :] = [cos(els.angle[i]) -sin(els.angle[i]) ; sin(els.angle[i]) cos(els.angle[i])]
        els.rotmatinv[i, :, :] = [cos(-els.angle[i]) -sin(-els.angle[i]) ; sin(-els.angle[i]) cos(-els.angle[i])]
        els.xnormal[i] = dy / magnitude
        els.ynormal[i] = -dx / magnitude
        els.xnodes[i, :] = [els.xcenter[i] - (2 / 3 * dx / 2), els.xcenter[i], els.xcenter[i] + (2 / 3 * dx / 2)]
        els.ynodes[i, :] = [els.ycenter[i] - (2 / 3 * dy / 2), els.ycenter[i], els.ycenter[i] + (2 / 3 * dy / 2)]
    end
    return nothing
end

export rotuσ
function rotuσ(u, σ, rotmat)
    # TODO: If this is slow hand expand the matrix vector μltiplies
    # inplace for speed.  Some benchmarks suggest 50x speedup!
    _u = zeros(size(u))
    _σ = zeros(size(σ))
    for i in 1:size(σ)[1]
        _u[i, 1], _u[i, 2] = rotmat * [u[i, 1] ; u[i, 2]]
        σtensor = [σ[i, 1] σ[i, 3] ; σ[i, 3] σ[i, 2]]
        _σ[i, 1], _σ[i, 3], _, _σ[i, 2] = rotmat * σtensor * rotmat'
    end
    return _u, _σ
end

export plotelements
function plotelements(els)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth=0.5)
    end
end

function stylesubplots(xlim, ylim)
    gca().set_aspect("equal")
    gca().set_xlim([xlim[1], xlim[2]])
    gca().set_ylim([ylim[1], ylim[2]])
    gca().set_xticks([xlim[1], xlim[2]])
    gca().set_yticks([ylim[1], ylim[2]])
end

function plotfields_contours(els, xobs, yobs, idx, field, title)
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
    plotelements(els)
end

export plotfields
function plotfields(els, xobs, yobs, disp, stress, suptitle)
    figure(figsize = (16, 8))
    subplot(2, 3, 1)
    quiver(xobs[:], yobs[:], disp[:, 1], disp[:, 2], units="width", color="b")
    stylesubplots([minimum(xobs) maximum(xobs)], [minimum(yobs) maximum(yobs)])
    plotelements(els)
    PyPlot.title("displacements")
    plotfields_contours(els, xobs, yobs, 2, disp[:, 1], L"u_x")
    plotfields_contours(els, xobs, yobs, 3, disp[:, 2], L"u_y")
    plotfields_contours(els, xobs, yobs, 4, stress[:, 1], L"\sigma_{xx}")
    plotfields_contours(els, xobs, yobs, 5, stress[:, 2], L"\sigma_{yy}")
    plotfields_contours(els, xobs, yobs, 6, stress[:, 3], L"\sigma_{xy}")
    PyPlot.suptitle(suptitle)
    tight_layout()
    show()
end

export ∂constslip
function ∂constslip(els, srcidx, obsidx, μ, nu)
    nobs = length(obsidx)
    nsrc = length(srcidx)
    ∂u = zeros(2 * nobs, 2 * nsrc)
    ∂σ = zeros(3 * nobs, 2 * nsrc)
    ∂t = zeros(2 * nobs, 2 * nsrc)

    for isrc in 1:nsrc
        for iobs in 1:nobs
            _∂u, _∂σ, _∂t = zeros(2, 2), zeros(3, 2), zeros(2, 2)
            _∂u[:, 1], _∂σ[:, 1] = constslip(
                els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]],
                els.halflength[srcidx[isrc]], μ, nu, 1, 0,
                els.xcenter[srcidx[isrc]], els.ycenter[srcidx[isrc]],
                els.rotmat[srcidx[isrc], :, :], els.rotmatinv[srcidx[isrc], :, :])
            _∂u[:, 2], _∂σ[:, 2] = constslip(
                els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]],
                els.halflength[srcidx[isrc]], μ, nu, 0, 1,
                els.xcenter[srcidx[isrc]], els.ycenter[srcidx[isrc]],
                els.rotmat[srcidx[isrc], :, :], els.rotmatinv[srcidx[isrc], :, :])
            _∂t[:, 1] = σ2t(_∂σ[:, 1], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            _∂t[:, 2] = σ2t(_∂σ[:, 2], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            ∂u[2*(iobs-1)+1 : 2*(iobs-1)+2, 2*(isrc-1)+1 : 2*(isrc-1)+2] = _∂u
            ∂σ[3*(iobs-1)+1 : 3*(iobs-1)+3, 2*(isrc-1)+1 : 2*(isrc-1)+2] = _∂σ
            ∂t[2*(iobs-1)+1 : 2*(iobs-1)+2, 2*(isrc-1)+1 : 2*(isrc-1)+2] = _∂t
        end
    end
    return ∂u, ∂σ, ∂t
end

export ∂consttrac
function ∂consttrac(els, srcidx, obsidx, μ, nu)
    nobs = length(obsidx)
    nsrc = length(srcidx)
    ∂u = zeros(2 * nobs, 2 * nsrc)
    ∂σ = zeros(3 * nobs, 2 * nsrc)
    ∂t = zeros(2 * nobs, 2 * nsrc)

    for isrc in 1:nsrc
        for iobs in 1:nobs
            _∂u, _∂σ, _∂t = zeros(2, 2), zeros(3, 2), zeros(2, 2)
            _∂u[:, 1], _∂σ[:, 1] = consttrac(
                els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]],
                els.halflength[srcidx[isrc]], μ, nu, 1, 0,
                els.xcenter[srcidx[isrc]], els.ycenter[srcidx[isrc]],
                els.rotmat[srcidx[isrc], :, :], els.rotmatinv[srcidx[isrc], :, :])
            _∂u[:, 2], _∂σ[:, 2] = consttrac(
                els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]],
                els.halflength[srcidx[isrc]], μ, nu, 0, 1,
                els.xcenter[srcidx[isrc]], els.ycenter[srcidx[isrc]],
                els.rotmat[srcidx[isrc], :, :], els.rotmatinv[srcidx[isrc], :, :])
            _∂t[:, 1] = σ2t(_∂σ[:, 1], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            _∂t[:, 2] = σ2t(_∂σ[:, 2], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            ∂u[2*(iobs-1)+1 : 2*(iobs-1)+2, 2*(isrc-1)+1 : 2*(isrc-1)+2] = _∂u
            ∂σ[3*(iobs-1)+1 : 3*(iobs-1)+3, 2*(isrc-1)+1 : 2*(isrc-1)+2] = _∂σ
            ∂t[2*(iobs-1)+1 : 2*(iobs-1)+2, 2*(isrc-1)+1 : 2*(isrc-1)+2] = _∂t
        end
    end
    return ∂u, ∂σ, ∂t
end

end
