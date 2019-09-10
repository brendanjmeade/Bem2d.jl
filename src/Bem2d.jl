module Bem2d
using LaTeXStrings
using PyCall
using PyPlot

export interleave
function interleave(vec1, vec2)
    return [vec1 vec2]'[:]
end

# Based on: http://julia-programming-language.2336112.n4.nabble.com/Meshgrid-function-td37003.html
export meshgrid
function meshgrid(xs, ys)
    return [xs[i] for i in 1:length(xs), j in 1:length(ys)], [ys[j] for i in 1:length(xs), j in 1:length(ys)]
end

# Create regular grid for plotting simple models
export obsgrid
function obsgrid(xmin, ymin, xmax, ymax, npts)
    xobs = LinRange(xmin, xmax, npts)
    yobs = LinRange(ymin, ymax, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]
    return xobs, yobs
end

export Elements
mutable struct Elements
    x1::Array{Float64,1}
    y1::Array{Float64,1}
    x2::Array{Float64,1}
    y2::Array{Float64,1}
    name::Array{String,1}
    uxconst::Array{Float64,1}
    uyconst::Array{Float64,1}
    uxquad::Array{Float64,2}
    uyquad::Array{Float64,2}
    angle::Array{Float64,1}
    length::Array{Float64,1}
    halflength::Array{Float64,1}
    xcenter::Array{Float64,1}
    ycenter::Array{Float64,1}
    rotmat::Array{Float64,3}
    rotmatinv::Array{Float64,3}
    xnormal::Array{Float64,1}
    ynormal::Array{Float64,1}
    xnodes::Array{Float64,2}
    ynodes::Array{Float64,2}
    a::Array{Float64,1}
    b::Array{Float64,1}
    σn::Array{Float64,1}
    endidx::Int64
    Elements(maxidx) = new(fill(NaN, maxidx), # x1
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
    x = range(xstart, stop = xend, length = npts)
    y = range(ystart, stop = yend, length = npts)
    x1 = x[1:1:end - 1]
    y1 = y[1:1:end - 1]
    x2 = x[2:1:end]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end

export multmatvec
function multmatvec(mats, vec1, vec2)
    newvec1 = zeros(size(vec1))
    newvec2 = zeros(size(vec2))
    for i in 1:length(vec1)
        newvec1[i], newvec2[i] = mats[i, :, :] * [vec1[i] ; vec2[i]]
    end
    return newvec1, newvec2
end

export multmatsinglevec
function multmatsinglevec(mat, vec1, vec2)
    newvec1 = zeros(size(vec1))
    newvec2 = zeros(size(vec2))
    for i in 1:length(vec1)
        newvec1[i], newvec2[i] = mat * [vec1[i] ; vec2[i]]
    end
    return newvec1, newvec2
end

export multmatvecsingle
function multmatvecsingle(mats, vec1, vec2)
    newvec1 = zeros(size(mats)[1])
    newvec2 = zeros(size(mats)[1])
    for i in 1:size(mats)[1]
        newvec1[i], newvec2[i] = mats[i, :, :] * [vec1 ; vec2]
    end
    return newvec1, newvec2
end

# Calculate u and σ for constant slip/traction elements
export constuσ
function constuσ(fun2uσ, x, y, els, idx, xcomp, ycomp, μ, ν)
    u, σ = zeros(length(x), 2), zeros(length(x), 3)
    for j in 1:length(idx)
        # Rotate and translate into SC coordinate system
        _x, _y = multmatsinglevec(els.rotmatinv[idx[j], :, :], x .- els.xcenter[idx[j]], y .- els.ycenter[idx[j]])
        _xcomp, _ycomp = els.rotmatinv[idx[j], :, :] * [xcomp[j] ; ycomp[j]]
        f = constkernel(_x, _y, els.halflength[idx[j]], ν)
        _u, _σ = fun2uσ(_xcomp, _ycomp, f, _y, μ, ν)
        _u, _σ = rotuσ(_u, _σ, els.rotmat[idx[j], :, :])
        u += _u
        σ += _σ
    end
    return u, σ
end

# Far-field displacements and stresses for constant quadratic elements
export quaduσ
function quaduσ(fun2uσ, x, y, els, idx, xcomp, ycomp, μ, ν)
    u, σ = zeros(length(x), 2), zeros(length(x), 3)
    for j in 1:length(idx)
        # Rotate and translate into SC coordinate system
        _x, _y = multmatsinglevec(els.rotmatinv[idx[j], :, :], x .- els.xcenter[idx[j]], y .- els.ycenter[idx[j]])
        f = quadkernel_farfield(_x, _y, els.halflength[idx[j]], ν)
        for i in 1:3
            _xcomp, _ycomp = els.rotmatinv[idx[j], :, :] * [xcomp[i] ; ycomp[i]]
            _u, _σ = fun2uσ(_xcomp, _ycomp, f[:, :, i], _y, μ, ν)
            _u, _σ = rotuσ(_u, _σ, els.rotmat[idx[j], :, :])
            u += _u
            σ += _σ
        end
    end
    return u, σ
end

# Far-field displacements and stresses for constant quadratic elements
export quaduσcoincident
function quaduσcoincident(fun2uσ, x, y, els, idx, xcomp, ycomp, μ, ν)
    u, σ = zeros(length(x), 2), zeros(length(x), 3)
    for j in 1:length(idx)
        # Rotate and translate into SC coordinate system
        _x, _y = multmatsinglevec(els.rotmatinv[idx[j], :, :], x .- els.xcenter[idx[j]], y .- els.ycenter[idx[j]])
        f = quadkernel_coincident(els.halflength[idx[j]], ν)
        for i in 1:3
            _xcomp, _ycomp = els.rotmatinv[idx[j], :, :] * [xcomp[i] ; ycomp[i]]
            _u, _σ = fun2uσ(_xcomp, _ycomp, f[:, :, i], _y, μ, ν)
            _u, _σ = rotuσ(_u, _σ, els.rotmat[idx[j], :, :])
            u += _u
            σ += _σ
        end
    end
    return u, σ
end


# Constant slip kernels from Starfield and Crouch, pages 49 and 82
function constkernel(x, y, a, ν)
    f = zeros(length(x), 7)
    for i in 1:length(x)
        f[i, 1] = -1 / (4 * π * (1 - ν)) * (y[i] * (atan(y[i], (x[i] - a)) - atan(y[i], (x[i] + a))) - (x[i] - a) * log(sqrt((x[i] - a)^2 + y[i]^2)) + (x[i] + a) * log(sqrt((x[i] + a)^2 + y[i]^2)))
        f[i, 2] = -1 / (4 * π * (1 - ν)) * ((atan(y[i], (x[i] - a)) - atan(y[i], (x[i] + a))))
        f[i, 3] = 1 / (4 * π * (1 - ν)) * (log(sqrt((x[i] - a)^2 + y[i]^2)) - log(sqrt((x[i] + a)^2 + y[i]^2)))
        f[i, 4] = 1 / (4 * π * (1 - ν)) * (y[i] / ((x[i] - a)^2 + y[i]^2) - y[i] / ((x[i] + a)^2 + y[i]^2))
        f[i, 5] = 1 / (4 * π * (1 - ν)) * ((x[i] - a) / ((x[i] - a)^2 + y[i]^2) - (x[i] + a) / ((x[i] + a)^2 + y[i]^2))
        f[i, 6] = 1 / (4 * π * (1 - ν)) * (((x[i] - a)^2 - y[i]^2) / ((x[i] - a)^2 + y[i]^2)^2 - ((x[i] + a)^2 - y[i]^2) / ((x[i] + a)^2 + y[i]^2)^2)
        f[i, 7] = 2 * y[i] / (4 * π * (1 - ν)) * ((x[i] - a) / ((x[i] - a)^2 + y[i]^2)^2 - (x[i] + a) / ((x[i] + a)^2 + y[i]^2)^2)
    end
    return f
end

# Coincident 3-node quadratic kernels
# Kernels for coincident integrals f, shape_function_idx, node_idx
export quadkernel_coincident
function quadkernel_coincident(a, ν)
    # TODO: Should I swap order: nodeidx, f, shapefunction Probably
    f = zeros(7, 3, 3)

    f[1, 1, 1] = -5 / 144 * a * log(25 / 9 * a^2) / (π - π * ν) - 17 / 288 * a * log(1 / 9 * a^2) / (π - π * ν) + 1 / 12 * a / (π - π * ν)
    f[1, 2, 1] = -25 / 288 * a * log(25 / 9 * a^2) / (π - π * ν) + 7 / 288 * a * log(1 / 9 * a^2) / (π - π * ν) + 1 / 12 * a / (π - π * ν)
    f[1, 3, 1] = -25 / 288 * a * log(25 / 9 * a^2) / (π - π * ν) - 1 / 144 * a * log(1 / 9 * a^2) / (π - π * ν) - 1 / 6 * a / (π - π * ν)
    f[1, 1, 2] = -3 / 16 * a * log(a) / (π - π * ν) - 1 / 8 * a / (π - π * ν)
    f[1, 2, 2] = -1 / 8 * a * log(a) / (π - π * ν) + 1 / 4 * a / (π - π * ν)
    f[1, 3, 2] = -3 / 16 * a * log(a) / (π - π * ν) - 1 / 8 * a / (π - π * ν)
    f[1, 1, 3] = -25 / 288 * a * log(25 / 9 * a^2) / (π - π * ν) - 1 / 144 * a * log(1 / 9 * a^2) / (π - π * ν) - 1 / 6 * a / (π - π * ν)
    f[1, 2, 3] = -25 / 288 * a * log(25 / 9 * a^2) / (π - π * ν) + 7 / 288 * a * log(1 / 9 * a^2) / (π - π * ν) + 1 / 12 * a / (π - π * ν)
    f[1, 3, 3] = -5 / 144 * a * log(25 / 9 * a^2) / (π - π * ν) - 17 / 288 * a * log(1 / 9 * a^2) / (π - π * ν) + 1 / 12 * a / (π - π * ν)

    f[2, 1, 1] = 1 / 4 / (ν - 1)
    f[2, 2, 1] = 0
    f[2, 3, 1] = 0
    f[2, 1, 2] = 0
    f[2, 2, 2] = 1 / 4 / (ν - 1)
    f[2, 3, 2] = 0
    f[2, 1, 3] = 0
    f[2, 2, 3] = 0
    f[2, 3, 3] = 1 / 4 / (ν - 1)

    f[3, 1, 1] = 1 / 8 * log(25 / 9 * a^2) / (π - π * ν) - 1 / 8 * log(1 / 9 * a^2) / (π - π * ν) - 3 / 4 / (π - π * ν)
    f[3, 2, 1] = 3 / 4 / (π - π * ν)
    f[3, 3, 1] = 0
    f[3, 1, 2] = -3 / 8 / (π - π * ν)
    f[3, 2, 2] = 0
    f[3, 3, 2] = 3 / 8 / (π - π * ν)
    f[3, 1, 3] = 0
    f[3, 2, 3] = -3 / 4 / (π - π * ν)
    f[3, 3, 3] = -1 / 8 * log(25 / 9 * a^2) / (π - π * ν) + 1 / 8 * log(1 / 9 * a^2) / (π - π * ν) + 3 / 4 / (π - π * ν)

    f[4, 1, 1] = -9 / 16 / (a * ν - a)
    f[4, 2, 1] = 3 / 4 / (a * ν - a)
    f[4, 3, 1] = -3 / 16 / (a * ν - a)
    f[4, 1, 2] = -3 / 16 / (a * ν - a)
    f[4, 2, 2] = 0
    f[4, 3, 2] = 3 / 16 / (a * ν - a)
    f[4, 1, 3] = 3 / 16 / (a * ν - a)
    f[4, 2, 3] = -3 / 4 / (a * ν - a)
    f[4, 3, 3] = 9 / 16 / (a * ν - a)

    f[5, 1, 1] = 9 / 32 * log(25 / 9 * a^2) / (π * a * ν - π * a) - 9 / 32 * log(1 / 9 * a^2) / (π * a * ν - π * a) + 27 / 80 / (π * a * ν - π * a)
    f[5, 2, 1] = -3 / 8 * log(25 / 9 * a^2) / (π * a * ν - π * a) + 3 / 8 * log(1 / 9 * a^2) / (π * a * ν - π * a) + 9 / 8 / (π * a * ν - π * a)
    f[5, 3, 1] = 3 / 32 * log(25 / 9 * a^2) / (π * a * ν - π * a) - 3 / 32 * log(1 / 9 * a^2) / (π * a * ν - π * a) - 9 / 16 / (π * a * ν - π * a)
    f[5, 1, 2] = -9 / 16 / (π * a * ν - π * a)
    f[5, 2, 2] = 13 / 8 / (π * a * ν - π * a)
    f[5, 3, 2] = -9 / 16 / (π * a * ν - π * a)
    f[5, 1, 3] = 3 / 32 * log(25 / 9 * a^2) / (π * a * ν - π * a) - 3 / 32 * log(1 / 9 * a^2) / (π * a * ν - π * a) - 9 / 16 / (π * a * ν - π * a)
    f[5, 2, 3] = -3 / 8 * log(25 / 9 * a^2) / (π * a * ν - π * a) + 3 / 8 * log(1 / 9 * a^2) / (π * a * ν - π * a) + 9 / 8 / (π * a * ν - π * a)
    f[5, 3, 3] = 9 / 32 * log(25 / 9 * a^2) / (π * a * ν - π * a) - 9 / 32 * log(1 / 9 * a^2) / (π * a * ν - π * a) + 27 / 80 / (π * a * ν - π * a)

    f[6, 1, 1] = 9 / 32 * log(25 / 9 * a^2) / (π * a^2 * ν - π * a^2) - 9 / 32 * log(1 / 9 * a^2) / (π * a^2 * ν - π * a^2) + 621 / 100 / (π * a^2 * ν - π * a^2)
    f[6, 2, 1] = -9 / 16 * log(25 / 9 * a^2) / (π * a^2 * ν - π * a^2) + 9 / 16 * log(1 / 9 * a^2) / (π * a^2 * ν - π * a^2) - 27 / 5 / (π * a^2 * ν - π * a^2)
    f[6, 3, 1] = 9 / 32 * log(25 / 9 * a^2) / (π * a^2 * ν - π * a^2) - 9 / 32 * log(1 / 9 * a^2) / (π * a^2 * ν - π * a^2) + 27 / 20 / (π * a^2 * ν - π * a^2)
    f[6, 1, 2] = 3 / 4 / (π * a^2 * ν - π * a^2)
    f[6, 2, 2] = 0
    f[6, 3, 2] = -3 / 4 / (π * a^2 * ν - π * a^2)
    f[6, 1, 3] = -9 / 32 * log(25 / 9 * a^2) / (π * a^2 * ν - π * a^2) + 9 / 32 * log(1 / 9 * a^2) / (π * a^2 * ν - π * a^2) - 27 / 20 / (π * a^2 * ν - π * a^2)
    f[6, 2, 3] = 9 / 16 * log(25 / 9 * a^2) / (π * a^2 * ν - π * a^2) - 9 / 16 * log(1 / 9 * a^2) / (π * a^2 * ν - π * a^2) + 27 / 5 / (π * a^2 * ν - π * a^2)
    f[6, 3, 3] = -9 / 32 * log(25 / 9 * a^2) / (π * a^2 * ν - π * a^2) + 9 / 32 * log(1 / 9 * a^2) / (π * a^2 * ν - π * a^2) - 621 / 100 / (π * a^2 * ν - π * a^2)

    f[7, 1, 1] = -9 / 16 / (a^2 * ν - a^2)
    f[7, 2, 1] = 9 / 8 / (a^2 * ν - a^2)
    f[7, 3, 1] = -9 / 16 / (a^2 * ν - a^2)
    f[7, 1, 2] = -9 / 16 / (a^2 * ν - a^2)
    f[7, 2, 2] = 9 / 8 / (a^2 * ν - a^2)
    f[7, 3, 2] = -9 / 16 / (a^2 * ν - a^2)
    f[7, 1, 3] = -9 / 16 / (a^2 * ν - a^2)
    f[7, 2, 3] = 9 / 8 / (a^2 * ν - a^2)
    f[7, 3, 3] = -9 / 16 / (a^2 * ν - a^2)

    f = permutedims(f, [3, 1, 2]) # This is just a hack until I get around to putting in the correct indices above
    return f
end

# Farfield 3-node quadratic kernels
# f has dimensions [nobs, f=7, shapefunctions=3]
export quadkernel_farfield
function quadkernel_farfield(x, y, a, ν)
    f = zeros(length(x), 7, 3)
    for i in 1:length(x)
        # term1 = atan((a - x[i]) / y[i])
        # term2 = atan((a + x[i]) / y[i])
        term1 = π / 2 * sign(y[i] / (a - x[i])) - atan(y[i] / (a - x[i]))
        term2 = π / 2 * sign(y[i] / (a + x[i])) - atan(y[i] / (a + x[i]))
        
        f[i, 1, 1] = -1 / 64 * ( 6 * y[i]^3 * (term2 + term1) - 6 * a^3 * log(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2) - 8 * a^3 - 12 * a^2 * x[i] + 12 * a * x[i]^2 - 12 * a * y[i]^2 + 6 * ( (2 * a * x[i] - 3 * x[i]^2) * term2 + (2 * a * x[i] - 3 * x[i]^2) * term1 ) * y[i] + 3 * (a * x[i]^2 - x[i]^3 - (a - 3 * x[i]) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 3 * (a * x[i]^2 - x[i]^3 - (a - 3 * x[i]) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) ) / (π * a^2 * ν - π * a^2)
        f[i, 1, 2] = 1 / 32 * (6 * y[i]^3 * (term2 + term1) + a^3 * log(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2) + a^3 * log(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2) - 8 * a^3 + 12 * a * x[i]^2 - 12 * a * y[i]^2 + 2 * ((4 * a^2 - 9 * x[i]^2) * term2 + (4 * a^2 - 9 * x[i]^2) * term1) * y[i] + (4 * a^2 * x[i] - 3 * x[i]^3 + 9 * x[i] * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (4 * a^2 * x[i] - 3 * x[i]^3 + 9 * x[i] * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^2 * ν - π * a^2)
        f[i, 1, 3] = -1 / 64 * (6 * y[i]^3 * (term2 + term1) - 6 * a^3 * log(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2) - 8 * a^3 + 12 * a^2 * x[i] + 12 * a * x[i]^2 - 12 * a * y[i]^2 - 6 * ((2 * a * x[i] + 3 * x[i]^2) * term2 + (2 * a * x[i] + 3 * x[i]^2) * term1) * y[i] - 3 * (a * x[i]^2 + x[i]^3 - (a + 3 * x[i]) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + 3 * (a * x[i]^2 + x[i]^3 - (a + 3 * x[i]) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^2 * ν - π * a^2)

        f[i, 2, 1] = -3 / 32 * (3 * y[i]^2 * (term2 + term1) - (a - 3 * x[i]) * y[i] * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (a - 3 * x[i]) * y[i] * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) - 6 * a * y[i] + (2 * a * x[i] - 3 * x[i]^2) * term2 + (2 * a * x[i] - 3 * x[i]^2) * term1) / (π * a^2 * ν - π * a^2)
        f[i, 2, 2] = 1 / 16 * (9 * y[i]^2 * (term2 + term1) + 9 * x[i] * y[i] * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * x[i] * y[i] * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) - 18 * a * y[i] + (4 * a^2 - 9 * x[i]^2) * term2 + (4 * a^2 - 9 * x[i]^2) * term1) / (π * a^2 * ν - π * a^2)
        f[i, 2, 3] = -3 / 32 * (3 * y[i]^2 * (term2 + term1) + (a + 3 * x[i]) * y[i] * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (a + 3 * x[i]) * y[i] * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) - 6 * a * y[i] - (2 * a * x[i] + 3 * x[i]^2) * term2 - (2 * a * x[i] + 3 * x[i]^2) * term1) / (π * a^2 * ν - π * a^2)

        f[i, 3, 1] = 3 / 64 * (8 * a^2 - 12 * a * x[i] - 4 * ((a - 3 * x[i]) * term2 + (a - 3 * x[i]) * term1) * y[i] - (2 * a * x[i] - 3 * x[i]^2 + 3 * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (2 * a * x[i] - 3 * x[i]^2 + 3 * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^2 * ν - π * a^2)
        f[i, 3, 2] = 1 / 32 * (36 * a * x[i] - 36 * (x[i] * term2 + x[i] * term1) * y[i] + (4 * a^2 - 9 * x[i]^2 + 9 * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (4 * a^2 - 9 * x[i]^2 + 9 * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^2 * ν - π * a^2)
        f[i, 3, 3] = -3 / 64 * (8 * a^2 + 12 * a * x[i] - 4 * ((a + 3 * x[i]) * term2 + (a + 3 * x[i]) * term1) * y[i] - (2 * a * x[i] + 3 * x[i]^2 - 3 * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (2 * a * x[i] + 3 * x[i]^2 - 3 * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^2 * ν - π * a^2)

        f[i, 4, 1] = 3 / 32 * (4 * a^2 * y[i]^3 - 2 * ((a - 3 * x[i]) * term2 + (a - 3 * x[i]) * term1) * y[i]^4 - 4 * ((a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * term2 + (a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * term1) * y[i]^2 + 4 * (a^4 - 3 * a^3 * x[i] + a^2 * x[i]^2) * y[i] - 2 * (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5) * term2 - 2 * (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5) * term1 - 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^6 * ν - π * a^6 + (π * a^2 * ν - π * a^2) * x[i]^4 + (π * a^2 * ν - π * a^2) * y[i]^4 - 2 * (π * a^4 * ν - π * a^4) * x[i]^2 + 2 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^2)
        f[i, 4, 2] = 1 / 16 * (20 * a^3 * x[i] * y[i] - 18 * (x[i] * term2 + x[i] * term1) * y[i]^4 - 36 * ((a^2 * x[i] + x[i]^3) * term2 + (a^2 * x[i] + x[i]^3) * term1) * y[i]^2 - 18 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5) * term2 - 18 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5) * term1 + 9 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^6 * ν - π * a^6 + (π * a^2 * ν - π * a^2) * x[i]^4 + (π * a^2 * ν - π * a^2) * y[i]^4 - 2 * (π * a^4 * ν - π * a^4) * x[i]^2 + 2 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^2)
        f[i, 4, 3] = -3 / 32 * (4 * a^2 * y[i]^3 - 2 * ((a + 3 * x[i]) * term2 + (a + 3 * x[i]) * term1) * y[i]^4 - 4 * ((a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * term2 + (a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * term1) * y[i]^2 + 4 * (a^4 + 3 * a^3 * x[i] + a^2 * x[i]^2) * y[i] - 2 * (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5) * term2 - 2 * (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5) * term1 + 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^6 * ν - π * a^6 + (π * a^2 * ν - π * a^2) * x[i]^4 + (π * a^2 * ν - π * a^2) * y[i]^4 - 2 * (π * a^4 * ν - π * a^4) * x[i]^2 + 2 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^2)

        f[i, 5, 1] = 3 / 32 * (6 * y[i]^5 * (term2 + term1) - 6 * a^5 - 4 * a^4 * x[i] + 18 * a^3 * x[i]^2 + 4 * a^2 * x[i]^3 - 12 * a * x[i]^4 - 12 * a * y[i]^4 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^3 - 2 * (9 * a^3 - 2 * a^2 * x[i] + 12 * a * x[i]^2) * y[i]^2 + 6 * ((a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term2 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term1) * y[i] - (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5 + (a - 3 * x[i]) * y[i]^4 + 2 * (a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5 + (a - 3 * x[i]) * y[i]^4 + 2 * (a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^6 * ν - π * a^6 + (π * a^2 * ν - π * a^2) * x[i]^4 + (π * a^2 * ν - π * a^2) * y[i]^4 - 2 * (π * a^4 * ν - π * a^4) * x[i]^2 + 2 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^2)
        f[i, 5, 2] = -1 / 16 * (18 * y[i]^5 * (term2 + term1) - 26 * a^5 + 62 * a^3 * x[i]^2 - 36 * a * x[i]^4 - 36 * a * y[i]^4 + 36 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^3 - 2 * (31 * a^3 + 36 * a * x[i]^2) * y[i]^2 + 18 * ((a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term2 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term1) * y[i] + 9 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5 + x[i] * y[i]^4 + 2 * (a^2 * x[i] + x[i]^3) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5 + x[i] * y[i]^4 + 2 * (a^2 * x[i] + x[i]^3) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^6 * ν - π * a^6 + (π * a^2 * ν - π * a^2) * x[i]^4 + (π * a^2 * ν - π * a^2) * y[i]^4 - 2 * (π * a^4 * ν - π * a^4) * x[i]^2 + 2 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^2)
        f[i, 5, 3] = 3 / 32 * (6 * y[i]^5 * (term2 + term1) - 6 * a^5 + 4 * a^4 * x[i] + 18 * a^3 * x[i]^2 - 4 * a^2 * x[i]^3 - 12 * a * x[i]^4 - 12 * a * y[i]^4 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^3 - 2 * (9 * a^3 + 2 * a^2 * x[i] + 12 * a * x[i]^2) * y[i]^2 + 6 * ((a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term2 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term1) * y[i] + (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5 + (a + 3 * x[i]) * y[i]^4 + 2 * (a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5 + (a + 3 * x[i]) * y[i]^4 + 2 * (a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^6 * ν - π * a^6 + (π * a^2 * ν - π * a^2) * x[i]^4 + (π * a^2 * ν - π * a^2) * y[i]^4 - 2 * (π * a^4 * ν - π * a^4) * x[i]^2 + 2 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^2)

        f[i, 6, 1] = 3 / 32 * (8 * a^8 - 24 * a^7 * x[i] - 16 * a^6 * x[i]^2 + 60 * a^5 * x[i]^3 + 8 * a^4 * x[i]^4 - 48 * a^3 * x[i]^5 + 12 * a * x[i]^7 + 12 * a * x[i] * y[i]^6 + 4 * (2 * a^4 + 12 * a^3 * x[i] + 9 * a * x[i]^3) * y[i]^4 + 4 * (4 * a^6 + 3 * a^5 * x[i] - 12 * a^4 * x[i]^2 + 9 * a * x[i]^5) * y[i]^2 - 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^10 * ν - π * a^10 + (π * a^2 * ν - π * a^2) * x[i]^8 + (π * a^2 * ν - π * a^2) * y[i]^8 - 4 * (π * a^4 * ν - π * a^4) * x[i]^6 + 4 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^6 + 6 * (π * a^6 * ν - π * a^6) * x[i]^4 + 2 * (3 * π * a^6 * ν - 3 * π * a^6 + 3 * (π * a^2 * ν - π * a^2) * x[i]^4 + 2 * (π * a^4 * ν - π * a^4) * x[i]^2) * y[i]^4 - 4 * (π * a^8 * ν - π * a^8) * x[i]^2 + 4 * (π * a^8 * ν - π * a^8 + (π * a^2 * ν - π * a^2) * x[i]^6 - (π * a^4 * ν - π * a^4) * x[i]^4 - (π * a^6 * ν - π * a^6) * x[i]^2) * y[i]^2)
        f[i, 6, 2] = 1 / 16 * (56 * a^7 * x[i] - 148 * a^5 * x[i]^3 + 128 * a^3 * x[i]^5 - 36 * a * x[i]^7 - 36 * a * x[i] * y[i]^6 - 12 * (8 * a^3 * x[i] + 9 * a * x[i]^3) * y[i]^4 - 4 * (a^5 * x[i] - 8 * a^3 * x[i]^3 + 27 * a * x[i]^5) * y[i]^2 + 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^10 * ν - π * a^10 + (π * a^2 * ν - π * a^2) * x[i]^8 + (π * a^2 * ν - π * a^2) * y[i]^8 - 4 * (π * a^4 * ν - π * a^4) * x[i]^6 + 4 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^6 + 6 * (π * a^6 * ν - π * a^6) * x[i]^4 + 2 * (3 * π * a^6 * ν - 3 * π * a^6 + 3 * (π * a^2 * ν - π * a^2) * x[i]^4 + 2 * (π * a^4 * ν - π * a^4) * x[i]^2) * y[i]^4 - 4 * (π * a^8 * ν - π * a^8) * x[i]^2 + 4 * (π * a^8 * ν - π * a^8 + (π * a^2 * ν - π * a^2) * x[i]^6 - (π * a^4 * ν - π * a^4) * x[i]^4 - (π * a^6 * ν - π * a^6) * x[i]^2) * y[i]^2)
        f[i, 6, 3] = -3 / 32 * (8 * a^8 + 24 * a^7 * x[i] - 16 * a^6 * x[i]^2 - 60 * a^5 * x[i]^3 + 8 * a^4 * x[i]^4 + 48 * a^3 * x[i]^5 - 12 * a * x[i]^7 - 12 * a * x[i] * y[i]^6 + 4 * (2 * a^4 - 12 * a^3 * x[i] - 9 * a * x[i]^3) * y[i]^4 + 4 * (4 * a^6 - 3 * a^5 * x[i] - 12 * a^4 * x[i]^2 - 9 * a * x[i]^5) * y[i]^2 + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (π * a^10 * ν - π * a^10 + (π * a^2 * ν - π * a^2) * x[i]^8 + (π * a^2 * ν - π * a^2) * y[i]^8 - 4 * (π * a^4 * ν - π * a^4) * x[i]^6 + 4 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^6 + 6 * (π * a^6 * ν - π * a^6) * x[i]^4 + 2 * (3 * π * a^6 * ν - 3 * π * a^6 + 3 * (π * a^2 * ν - π * a^2) * x[i]^4 + 2 * (π * a^4 * ν - π * a^4) * x[i]^2) * y[i]^4 - 4 * (π * a^8 * ν - π * a^8) * x[i]^2 + 4 * (π * a^8 * ν - π * a^8 + (π * a^2 * ν - π * a^2) * x[i]^6 - (π * a^4 * ν - π * a^4) * x[i]^4 - (π * a^6 * ν - π * a^6) * x[i]^2) * y[i]^2)

        f[i, 7, 1] = -3 / 16 * (3 * y[i]^8 * (term2 + term1) - 6 * a * y[i]^7 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^6 - 6 * (4 * a^3 + 3 * a * x[i]^2) * y[i]^5 + 6 * ((3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term2 + (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term1) * y[i]^4 - 2 * (15 * a^5 - 8 * a^4 * x[i] + 9 * a * x[i]^4) * y[i]^3 + 12 * ((a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term2 + (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term1) * y[i]^2 - 2 * (6 * a^7 - 8 * a^6 * x[i] + 3 * a^5 * x[i]^2 + 8 * a^4 * x[i]^3 - 12 * a^3 * x[i]^4 + 3 * a * x[i]^6) * y[i] + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term2 + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term1) / (π * a^10 * ν - π * a^10 + (π * a^2 * ν - π * a^2) * x[i]^8 + (π * a^2 * ν - π * a^2) * y[i]^8 - 4 * (π * a^4 * ν - π * a^4) * x[i]^6 + 4 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^6 + 6 * (π * a^6 * ν - π * a^6) * x[i]^4 + 2 * (3 * π * a^6 * ν - 3 * π * a^6 + 3 * (π * a^2 * ν - π * a^2) * x[i]^4 + 2 * (π * a^4 * ν - π * a^4) * x[i]^2) * y[i]^4 - 4 * (π * a^8 * ν - π * a^8) * x[i]^2 + 4 * (π * a^8 * ν - π * a^8 + (π * a^2 * ν - π * a^2) * x[i]^6 - (π * a^4 * ν - π * a^4) * x[i]^4 - (π * a^6 * ν - π * a^6) * x[i]^2) * y[i]^2)
        f[i, 7, 2] = 1 / 8 * (9 * y[i]^8 * (term2 + term1) - 18 * a * y[i]^7 + 36 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^6 - 2 * (32 * a^3 + 27 * a * x[i]^2) * y[i]^5 + 18 * ((3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term2 + (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term1) * y[i]^4 - 2 * (37 * a^5 + 8 * a^3 * x[i]^2 + 27 * a * x[i]^4) * y[i]^3 + 36 * ((a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term2 + (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term1) * y[i]^2 - 2 * (14 * a^7 + a^5 * x[i]^2 - 24 * a^3 * x[i]^4 + 9 * a * x[i]^6) * y[i] + 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term2 + 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term1) / (π * a^10 * ν - π * a^10 + (π * a^2 * ν - π * a^2) * x[i]^8 + (π * a^2 * ν - π * a^2) * y[i]^8 - 4 * (π * a^4 * ν - π * a^4) * x[i]^6 + 4 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^6 + 6 * (π * a^6 * ν - π * a^6) * x[i]^4 + 2 * (3 * π * a^6 * ν - 3 * π * a^6 + 3 * (π * a^2 * ν - π * a^2) * x[i]^4 + 2 * (π * a^4 * ν - π * a^4) * x[i]^2) * y[i]^4 - 4 * (π * a^8 * ν - π * a^8) * x[i]^2 + 4 * (π * a^8 * ν - π * a^8 + (π * a^2 * ν - π * a^2) * x[i]^6 - (π * a^4 * ν - π * a^4) * x[i]^4 - (π * a^6 * ν - π * a^6) * x[i]^2) * y[i]^2)
        f[i, 7, 3] = -3 / 16 * (3 * y[i]^8 * (term2 + term1) - 6 * a * y[i]^7 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^6 - 6 * (4 * a^3 + 3 * a * x[i]^2) * y[i]^5 + 6 * ((3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term2 + (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term1) * y[i]^4 - 2 * (15 * a^5 + 8 * a^4 * x[i] + 9 * a * x[i]^4) * y[i]^3 + 12 * ((a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term2 + (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term1) * y[i]^2 - 2 * (6 * a^7 + 8 * a^6 * x[i] + 3 * a^5 * x[i]^2 - 8 * a^4 * x[i]^3 - 12 * a^3 * x[i]^4 + 3 * a * x[i]^6) * y[i] + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term2 + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term1) / (π * a^10 * ν - π * a^10 + (π * a^2 * ν - π * a^2) * x[i]^8 + (π * a^2 * ν - π * a^2) * y[i]^8 - 4 * (π * a^4 * ν - π * a^4) * x[i]^6 + 4 * (π * a^4 * ν - π * a^4 + (π * a^2 * ν - π * a^2) * x[i]^2) * y[i]^6 + 6 * (π * a^6 * ν - π * a^6) * x[i]^4 + 2 * (3 * π * a^6 * ν - 3 * π * a^6 + 3 * (π * a^2 * ν - π * a^2) * x[i]^4 + 2 * (π * a^4 * ν - π * a^4) * x[i]^2) * y[i]^4 - 4 * (π * a^8 * ν - π * a^8) * x[i]^2 + 4 * (π * a^8 * ν - π * a^8 + (π * a^2 * ν - π * a^2) * x[i]^6 - (π * a^4 * ν - π * a^4) * x[i]^4 - (π * a^6 * ν - π * a^6) * x[i]^2) * y[i]^2)
    end
    return f
end

# Generalization from Starfield and Crouch
export trac2uσ
function trac2uσ(xcomp, ycomp, f, y, μ, ν)
    u, σ = zeros(length(y), 2), zeros(length(y), 3)
    _xcomp, _ycomp = -xcomp, -ycomp # For Okada consistency
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
    u, σ = zeros(length(y), 2), zeros(length(y), 3)
    _xcomp, _ycomp = -xcomp, -ycomp # For Okada consistency
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
    _u, _σ = zeros(size(u)), zeros(size(σ))
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
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 0.5)
    end
    return nothing
end

function stylesubplots(xlim, ylim)
    gca().set_aspect("equal")
    gca().set_xlim([xlim[1], xlim[2]])
    gca().set_ylim([ylim[1], ylim[2]])
    gca().set_xticks([xlim[1], xlim[2]])
    gca().set_yticks([ylim[1], ylim[2]])
    return nothing
end

function plotfields_contours(els, xobs, yobs, idx, field, title)
    ncontours = 20
    xlim = [minimum(xobs) maximum(xobs)]
    ylim = [minimum(yobs) maximum(yobs)]
    subplot(2, 3, idx)
    scale = 5e-1
    fieldmax = maximum(@.abs(field))
    contourf(xobs, yobs, reshape(field, size(xobs)), ncontours,
        vmin = -scale * fieldmax, vmax = scale * fieldmax, cmap = plt.get_cmap("PiYG"))
    clim(-scale * fieldmax, scale * fieldmax)
    colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    contour(xobs, yobs, reshape(field, size(xobs)), ncontours,
        vmin = -scale * fieldmax, vmax = scale * fieldmax, linewidths = 0.25, colors = "k")
    PyPlot.title(title)
    stylesubplots(xlim, ylim)
    plotelements(els)
    return nothing
end

export plotfields
function plotfields(els, xobs, yobs, disp, stress, suptitle)
    figure(figsize = (16, 8))
    subplot(2, 3, 1)
    quiver(xobs[:], yobs[:], disp[:, 1], disp[:, 2], units = "width", color = "b")
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
    return nothing
end

export ∂constuσ
function ∂constuσ(fun2uσ, els, srcidx, obsidx, μ, ν)
    nobs, nsrc = length(obsidx), length(srcidx)
    ∂u, ∂σ, ∂t = zeros(2 * nobs, 2 * nsrc), zeros(3 * nobs, 2 * nsrc), zeros(2 * nobs, 2 * nsrc)
    for isrc in 1:nsrc
        for iobs in 1:nobs
            _∂u, _∂σ, _∂t = zeros(2, 2), zeros(3, 2), zeros(2, 2)
            _∂u[:, 1], _∂σ[:, 1] = constuσ(fun2uσ, els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]], els, srcidx[isrc], 1, 0, μ, ν)
            _∂u[:, 2], _∂σ[:, 2] = constuσ(fun2uσ, els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]], els, srcidx[isrc], 0, 1, μ, ν)
            _∂t[:, 1] = σ2t(_∂σ[:, 1], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            _∂t[:, 2] = σ2t(_∂σ[:, 2], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            ∂u[2 * (iobs - 1) + 1:2 * (iobs - 1) + 2, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _∂u
            ∂σ[3 * (iobs - 1) + 1:3 * (iobs - 1) + 3, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _∂σ
            ∂t[2 * (iobs - 1) + 1:2 * (iobs - 1) + 2, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _∂t
        end
    end
    return ∂u, ∂σ, ∂t
end

# This is small helper function to return interleaved displacements
# and stresses rather than the stacked columns
export quaduσinterleaved
function quaduσinterleaved(fun2uσ, xobs, yobs, els, idx, para, perp, μ, ν)
    ustacked, σstacked = quaduσ(fun2uσ, xobs, yobs, els, idx, para, perp, μ, ν)
    uinterleaved = zeros(2 * length(xobs))
    uinterleaved[1:2:end] = ustacked[:, 1]
    uinterleaved[2:2:end] = ustacked[:, 2]
    σinterleaved = zeros(3 * length(xobs))
    σinterleaved[1:3:end] = σstacked[:, 1]
    σinterleaved[2:3:end] = σstacked[:, 2]
    σinterleaved[3:3:end] = σstacked[:, 3]
    return uinterleaved, σinterleaved
end

export quaduσcoincidentinterleaved
function quaduσcoincidentinterleaved(fun2uσ, xobs, yobs, els, idx, para, perp, μ, ν)
    ustacked, σstacked = quaduσcoincident(fun2uσ, xobs, yobs, els, idx, para, perp, μ, ν)
    uinterleaved = zeros(2 * length(xobs))
    uinterleaved[1:2:end] = ustacked[:, 1]
    uinterleaved[2:2:end] = ustacked[:, 2]
    σinterleaved = zeros(3 * length(xobs))
    σinterleaved[1:3:end] = σstacked[:, 1]
    σinterleaved[2:3:end] = σstacked[:, 2]
    σinterleaved[3:3:end] = σstacked[:, 3]
    return uinterleaved, σinterleaved
end

# TODO.  This is under active development
# Each 6x6 displacement submatrix is shaped as:
#
#       sx(ϕ1) sy(ϕ1) sx(ϕ2) sy(ϕ2) sx(ϕ3) sy(ϕ3)
# ux(ϕ1)
# uy(ϕ1)
# ux(ϕ2)
# uy(ϕ2)
# ux(ϕ3)
# uy(ϕ3)
export ∂quaduσ
function ∂quaduσ(fun2uσ, els, srcidx, obsidx, μ, ν)
    nobs, nsrc = length(obsidx), length(srcidx)
    ∂u, ∂σ, ∂t = zeros(6 * nobs, 6 * nsrc), zeros(9 * nobs, 6 * nsrc), zeros(6 * nobs, 6 * nsrc)

    # Far-field interactions for all except main diagonal
    for isrc in 1:nsrc
        for iobs in 1:nobs
            _∂u, _∂σ, _∂t = zeros(6, 6), zeros(9, 6), zeros(6, 6)
            if isrc == iobs
                # println("coincident integrals")
                _∂u[:, 1], _∂σ[:, 1] = quaduσcoincidentinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [1 0 0], [0 0 0], μ, ν)
                _∂u[:, 2], _∂σ[:, 2] = quaduσcoincidentinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [1 0 0], μ, ν)
                _∂u[:, 3], _∂σ[:, 3] = quaduσcoincidentinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 1 0], [0 0 0], μ, ν)
                _∂u[:, 4], _∂σ[:, 4] = quaduσcoincidentinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 1 0], μ, ν)
                _∂u[:, 5], _∂σ[:, 5] = quaduσcoincidentinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 1], [0 0 0], μ, ν)
                _∂u[:, 6], _∂σ[:, 6] = quaduσcoincidentinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 0 1], μ, ν)
            else
                # println("farfield integrals")
                _∂u[:, 1], _∂σ[:, 1] = quaduσinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [1 0 0], [0 0 0], μ, ν)
                _∂u[:, 2], _∂σ[:, 2] = quaduσinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [1 0 0], μ, ν)
                _∂u[:, 3], _∂σ[:, 3] = quaduσinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 1 0], [0 0 0], μ, ν)
                _∂u[:, 4], _∂σ[:, 4] = quaduσinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 1 0], μ, ν)
                _∂u[:, 5], _∂σ[:, 5] = quaduσinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 1], [0 0 0], μ, ν)
                _∂u[:, 6], _∂σ[:, 6] = quaduσinterleaved(fun2uσ, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 0 1], μ, ν)
            end
            # _∂t[:, 1] = σ2t(_∂σ[:, 1], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            # _∂t[:, 2] = σ2t(_∂σ[:, 2], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            ∂u[6 * (iobs - 1) + 1:6 * (iobs - 1) + 6, 6 * (isrc - 1) + 1:6 * (isrc - 1) + 6] = _∂u
            ∂σ[9 * (iobs - 1) + 1:9 * (iobs - 1) + 9, 6 * (isrc - 1) + 1:6 * (isrc - 1) + 6] = _∂σ
            ∂t[6 * (iobs - 1) + 1:6 * (iobs - 1) + 6, 6 * (isrc - 1) + 1:6 * (isrc - 1) + 6] = _∂t
        end
    end

    # Coincident interactions for self interactions on main diagonal

    return ∂u, ∂σ, ∂t
end


end
x₁₍