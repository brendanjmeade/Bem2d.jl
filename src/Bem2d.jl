module Bem2d
using LaTeXStrings
using Colors
using ColorSchemes
using PyCall
using PyPlot
using Infiltrator


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
    normalstress::Array{Float64,1}
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


export discretizedarc
"""
    discretized_arc(thetastart, thetaend, radius, n_pts)

Generate regularly spaced eleemnts along an curved arc.
"""
function discretizedarc(thetastart, thetaend, radius, n_pts)
    # Create geometry of discretized arc
    thetarange = collect(LinRange(thetastart, thetaend, n_pts + 1))
    x = @. radius * cos(thetarange)
    y = @. radius * sin(thetarange)
    x1 = x[1:1:end-1]
    x2 = x[2:1:end]
    y1 = y[1:1:end-1]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end


export addelsez!
"""
    addelsez!()

Covenience function for quicly adding elements with geometry and name only
"""
function addelsez!(els, x1, y1, x2, y2, namestring)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = namestring
    end
    standardize_elements!(els)
end


export flipud
"""
    flipud()

Just an alias to the 1-linear that replicates matlab's matrix flipud
Taken from:
https://cheatsheets.quantecon.org/#manipulating-vectors-and-matrices
"""
function flipud(mat)
    mat = reverse(mat, dims = 1)
    return(mat)
end


export interleave
"""
    interleave(vec1, vec2)

Interleave values in two arrays of equal length.

# Example
```julia-repl
julia> interleave([1, 2], [3, 4])
4-element Array{Int64,1}:
 1
 3
 2
 4
```
"""
function interleave(vec1, vec2)
    return transpose([vec1 vec2])[:]
end


export meshgrid
"""
    meshgrid(xs, ys)

Replicate Matlab's meshgrid.  Makes plaid grids from two vectors.
"""
function meshgrid(xs, ys)
    xmat = repeat(xs, 1, length(ys))
    ymat = transpose(repeat(ys, 1, length(xs)))
    return xmat, ymat
end


export obsgrid
"""
    obsgrid(xmin, ymin, xmax, ymax, npts)

Create set of observation corrdinates on a regular rectangular grid.
"""
function obsgrid(xmin, ymin, xmax, ymax, npts)
    T = typeof(xmin)
    xobs, yobs = meshgrid(LinRange(T(xmin), T(xmax), npts), LinRange(T(ymin), T(ymax), npts))
    return xobs[:], yobs[:]
end


export updateendidx!
"""
    updateendidx!(els)

In place determinate of last index of structures els in which there is data
"""
function updateendidx!(els)
    els.endidx = findall(isnan, els.x1)[1] - 1
    return nothing
end


export discretizedline
"""
    discretizedline(xstart, ystart, xend, yend, nelements)

Calculated start and endpoint coordinates for a line discretized 
into nelements
"""
function discretizedline(xstart, ystart, xend, yend, nelements)
    npts = nelements + 1
    x = range(xstart, stop = xend, length = npts)
    y = range(ystart, stop = yend, length = npts)
    x1 = x[1:1:end - 1]
    y1 = y[1:1:end - 1]
    x2 = x[2:1:end]
    y2 = y[2:1:end]
    return collect(x1), collect(y1), collect(x2), collect(y2)
end


export multmatvec
"""
    multmatvec(mats, vec1, vec2)

Matrix vector multiply. Multiple matrices.  Multiple vectors. 
"""
function multmatvec(mats, vec1, vec2)
    newvec1 = zeros(size(vec1))
    newvec2 = zeros(size(vec2))
    for i in 1:length(vec1)
        newvec1[i], newvec2[i] = mats[i, :, :] * [vec1[i] ; vec2[i]]
    end
    return newvec1, newvec2
end


export multmatvec!
"""
    multmatvec!(mats, vec1, vec2)

In place matrix vector multiply. Multiple matrices.  Multiple vectors. 
"""
function multmatvec!(vx, vy, mats, vec1, vec2)
    @inbounds for i in 1:length(vec1)
        @views vx[i] = mats[i, 1, 1] * vec1[i] + mats[i, 1, 2] * vec2[i]
        @views vy[i] = mats[i, 2, 1] * vec1[i] + mats[i, 2, 2] * vec2[i]
    end
    return nothing
end


export multmatvecquad!
"""
    multmatvecquad!(mats, vec1, vec2)

In place matrix vector multiply. Multiple matrices.  Multiple vectors. 
Built for quadratic case.
"""
function multmatvecquad!(vx, vy, mats, vec1, vec2)
    for i in 1:length(vec1)
        matidx = Int64(floor((i - 1) / 3) + 1) # Change w/ every 3rd node
        @views vx[i] = mats[matidx, 1, 1] * vec1[i] + mats[matidx, 1, 2] * vec2[i]
        @views vy[i] = mats[matidx, 2, 1] * vec1[i] + mats[matidx, 2, 2] * vec2[i]
    end
    return nothing
end


export multmatsinglevec
"""
    multmatvecsingle(mats, vec1, vec2)

Matrix vector multiply. Single matrix.  Multiple vectors. 
"""
function multmatsinglevec(mat, vec1, vec2)
    newvec1 = zeros(size(vec1))
    newvec2 = zeros(size(vec2))
    for i in 1:length(vec1)
        newvec1[i], newvec2[i] = mat * [vec1[i] ; vec2[i]]
    end
    return newvec1, newvec2
end


export multmatvecsingle
"""
    multmatvecsingle(mats, vec1, vec2)

Matrix vector multiply. Multiple matrices.  Single set of vectors. 
"""
function multmatvecsingle(mats, vec1, vec2)
    newvec1 = zeros(size(mats)[1])
    newvec2 = zeros(size(mats)[1])
    for i in 1:size(mats)[1]
        newvec1[i], newvec2[i] = mats[i, :, :] * [vec1 ; vec2]
    end
    return newvec1, newvec2
end


export constdispstress
"""
    constdispstress(fun2dispstress, x, y, els, idx, xcomp, ycomp, mu, nu)

Calculate displacements in global reference frame for constant elements.
"""
function constdispstress(fun2dispstress, x, y, els, idx, xcomp, ycomp, mu, nu)
    disp, stress = zeros(length(x), 2), zeros(length(x), 3)
    _x, _y = zeros(length(x)), zeros(length(x))
    _xcomp, _ycomp = zeros(1), zeros(1)
    f = zeros(length(x), 7)
    _disp, _stress = zeros(length(x), 2), zeros(length(x), 3)
    @inbounds @simd for j in 1:length(idx)
        @views _x, _y = multmatsinglevec(els.rotmatinv[idx[j], :, :], x .- els.xcenter[idx[j]], y .- els.ycenter[idx[j]])
        @views _xcomp, _ycomp = els.rotmatinv[idx[j], :, :] * [xcomp[j] ; ycomp[j]]
        @views f = constkernel(_x, _y, els.halflength[idx[j]], nu)
        _disp, _stress = fun2dispstress(_xcomp, _ycomp, f, _y, mu, nu)
        @views _disp, _stress = rotdispstress(_disp, _stress, els.rotmat[idx[j], :, :])
        disp += _disp
        stress += _stress
    end
    return disp, stress
end


export quaddispstress
"""
    quaddispstress(fun2dispstress, x, y, els, idx, xcomp, ycomp, mu, nu)

Calculate displacements in global reference frame for quadratic elements.
Far-field case.
"""
function quaddispstress(fun2dispstress, x, y, els, idx, xcomp, ycomp, mu, nu)
    disp = zeros(length(x), 2)
    stress = zeros(length(x), 3)

    for j in 1:length(idx)
        _x, _y = multmatsinglevec(els.rotmatinv[idx[j], :, :], x .- els.xcenter[idx[j]], y .- els.ycenter[idx[j]])
        f = quadkernel_farfield(_x, _y, els.halflength[idx[j]], nu)
        _xcomp, _ycomp = multmatsinglevec(els.rotmatinv[idx[j], :, :], xcomp[j, :], ycomp[j, :])

        # TODO: I think I need to get xnodes in a local coordinate system
        # This should be a translate and a rotate
        # Should these matrices be moved to standardize_elements?
        xnodes, ynodes = multmatsinglevec(els.rotmatinv[idx[j], :, :], els.xnodes[idx[j], :] .- els.xcenter[idx[j]], els.ynodes[idx[j], :] .- els.ycenter[idx[j]])
        phix = slip2coef(xnodes, _xcomp, els.halflength[idx[j]])
        phiy = slip2coef(xnodes, _ycomp, els.halflength[idx[j]])
        for i in 1:3
            # _disp, _stress = fun2dispstress(_xcomp[i], _ycomp[i], f[:, :, i], _y, mu, nu)
            _disp, _stress = fun2dispstress(phix[i], phiy[i], f[:, :, i], _y, mu, nu)
            _disp, _stress = rotdispstress(_disp, _stress, els.rotmat[idx[j], :, :])
            disp += _disp
            stress += _stress
        end
    end
    return disp, stress
end


export quaddispstresscoincident
"""
    quaddispstresscoincident(fun2dispstress, x, y, els, idx, xcomp, ycomp, mu, nu)

Calculate displacements in global reference frame for quadratic elements.
Coincident case.
"""
function quaddispstresscoincident(fun2dispstress, x, y, els, idx, xcomp, ycomp, mu, nu)
    disp = zeros(length(x), 2)
    stress = zeros(length(x), 3)
    for j in 1:length(idx)
        # Rotate and translate into SC coordinate system
        _x, _y = multmatsinglevec(els.rotmatinv[idx[j], :, :], x .- els.xcenter[idx[j]], y .- els.ycenter[idx[j]])
        f = quadkernel_coincident(els.halflength[idx[j]], nu)
        _xcomp, _ycomp = multmatsinglevec(els.rotmatinv[idx[j], :, :], xcomp[j, :], ycomp[j, :])
        xnodes, ynodes = multmatsinglevec(els.rotmatinv[idx[j], :, :], els.xnodes[idx[j], :] .- els.xcenter[idx[j]], els.ynodes[idx[j], :] .- els.ycenter[idx[j]])
        phix = slip2coef(xnodes, _xcomp, els.halflength[idx[j]])
        phiy = slip2coef(xnodes, _ycomp, els.halflength[idx[j]])
        for i in 1:3
            # _disp, _stress = fun2uσ(_xcomp[i], _ycomp[i], f[:, :, i], _y, mu, nu)
            _disp, _stress = fun2dispstress(phix[i], phiy[i], f[:, :, i], _y, mu, nu)
            _disp, _stress = rotdispstress(_disp, _stress, els.rotmat[idx[j], :, :])
            disp += _disp
            stress += _stress
        end
    end
    return disp, stress
end


"""
    constkernel(x, y, a, nu)

Calculate constant slip kernels from Starfield and Crouch, pages 49 and 82.
"""
function constkernel(x, y, a, nu)
    f = zeros(length(x), 7)
    leadingconst = 1/(4*pi*(1-nu))
    for i in 1:length(x)
        # f(x) = f_1
        f[i, 1] = -leadingconst * (y[i] * (atan(y[i], (x[i]-a)) - atan(y[i], (x[i]+a))) - (x[i] - a) * log(sqrt((x[i] - a)^2 + y[i]^2)) + (x[i]+a) * log(sqrt((x[i]+a)^2 + y[i]^2)))

        # df/dy = f_2
        f[i, 2] = -leadingconst * ((atan(y[i], (x[i]-a)) - atan(y[i], (x[i]+a))))

        # df/dx = f_3
        f[i, 3] = leadingconst * (log(sqrt((x[i]-a)^2 + y[i]^2)) - log(sqrt((x[i]+a)^2 + y[i]^2)))

        # d2f/dxdy = f_4
        f[i, 4] = leadingconst * (y[i] / ((x[i]-a)^2 + y[i]^2) - y[i] / ((x[i]+a)^2 + y[i]^2))

        # d2f/dxdx = -d2f/dydy = f_5
        f[i, 5] = leadingconst * ((x[i]-a) / ((x[i]-a)^2 + y[i]^2) - (x[i]+a) / ((x[i]+a)^2 + y[i]^2))

        # d3f/dxdydy = -d3f/dxdxdx = f_6
        f[i, 6] = leadingconst * (((x[i]-a)^2 - y[i]^2) / ((x[i]-a)^2 + y[i]^2)^2 - ((x[i]+a)^2 - y[i]^2) / ((x[i]+a)^2 + y[i]^2)^2)

        # d3f/dydydy = -d3f/dxdxdy = f7
        f[i, 7] = 2 * y[i] / (4 * pi * (1 - nu)) * ((x[i]-a) / ((x[i]-a)^2 + y[i]^2)^2 - (x[i]+a) / ((x[i]+a)^2 + y[i]^2)^2)
    end
    return f
end


export quadkernel_coincident
"""
    quadkernel_coincident(a, nu)
    
Farfield 3-node quadratic kernels for substitution into CS kernels
f has dimensions [nobs, f=7, shapefunctions=3]

Note: These were calculated symbolically using Maxima and then
reduced to eliminate singularites by zeroing out diverging terms.
This was very tedious.
"""
function quadkernel_coincident(a, nu)
    # TODO: Should I swap order: nodeidx, f, shapefunction Probably
    f = zeros(7, 3, 3)

    f[1, 1, 1] = -5 / 144 * a * log(25 / 9 * a^2) / (pi - pi * nu) - 17 / 288 * a * log(1 / 9 * a^2) / (pi - pi * nu) + 1 / 12 * a / (pi - pi * nu)
    f[1, 2, 1] = -25 / 288 * a * log(25 / 9 * a^2) / (pi - pi * nu) + 7 / 288 * a * log(1 / 9 * a^2) / (pi - pi * nu) + 1 / 12 * a / (pi - pi * nu)
    f[1, 3, 1] = -25 / 288 * a * log(25 / 9 * a^2) / (pi - pi * nu) - 1 / 144 * a * log(1 / 9 * a^2) / (pi - pi * nu) - 1 / 6 * a / (pi - pi * nu)
    f[1, 1, 2] = -3 / 16 * a * log(a) / (pi - pi * nu) - 1 / 8 * a / (pi - pi * nu)
    f[1, 2, 2] = -1 / 8 * a * log(a) / (pi - pi * nu) + 1 / 4 * a / (pi - pi * nu)
    f[1, 3, 2] = -3 / 16 * a * log(a) / (pi - pi * nu) - 1 / 8 * a / (pi - pi * nu)
    f[1, 1, 3] = -25 / 288 * a * log(25 / 9 * a^2) / (pi - pi * nu) - 1 / 144 * a * log(1 / 9 * a^2) / (pi - pi * nu) - 1 / 6 * a / (pi - pi * nu)
    f[1, 2, 3] = -25 / 288 * a * log(25 / 9 * a^2) / (pi - pi * nu) + 7 / 288 * a * log(1 / 9 * a^2) / (pi - pi * nu) + 1 / 12 * a / (pi - pi * nu)
    f[1, 3, 3] = -5 / 144 * a * log(25 / 9 * a^2) / (pi - pi * nu) - 17 / 288 * a * log(1 / 9 * a^2) / (pi - pi * nu) + 1 / 12 * a / (pi - pi * nu)

    f[2, 1, 1] = 1 / 4 / (nu - 1)
    f[2, 2, 1] = 0
    f[2, 3, 1] = 0
    f[2, 1, 2] = 0
    f[2, 2, 2] = 1 / 4 / (nu - 1)
    f[2, 3, 2] = 0
    f[2, 1, 3] = 0
    f[2, 2, 3] = 0
    f[2, 3, 3] = 1 / 4 / (nu - 1)

    f[3, 1, 1] = 1 / 8 * log(25 / 9 * a^2) / (pi - pi * nu) - 1 / 8 * log(1 / 9 * a^2) / (pi - pi * nu) - 3 / 4 / (pi - pi * nu)
    f[3, 2, 1] = 3 / 4 / (pi - pi * nu)
    f[3, 3, 1] = 0
    f[3, 1, 2] = -3 / 8 / (pi - pi * nu)
    f[3, 2, 2] = 0
    f[3, 3, 2] = 3 / 8 / (pi - pi * nu)
    f[3, 1, 3] = 0
    f[3, 2, 3] = -3 / 4 / (pi - pi * nu)
    f[3, 3, 3] = -1 / 8 * log(25 / 9 * a^2) / (pi - pi * nu) + 1 / 8 * log(1 / 9 * a^2) / (pi - pi * nu) + 3 / 4 / (pi - pi * nu)

    f[4, 1, 1] = -9 / 16 / (a * nu - a)
    f[4, 2, 1] = 3 / 4 / (a * nu - a)
    f[4, 3, 1] = -3 / 16 / (a * nu - a)
    f[4, 1, 2] = -3 / 16 / (a * nu - a)
    f[4, 2, 2] = 0
    f[4, 3, 2] = 3 / 16 / (a * nu - a)
    f[4, 1, 3] = 3 / 16 / (a * nu - a)
    f[4, 2, 3] = -3 / 4 / (a * nu - a)
    f[4, 3, 3] = 9 / 16 / (a * nu - a)

    f[5, 1, 1] = 9 / 32 * log(25 / 9 * a^2) / (pi * a * nu - pi * a) - 9 / 32 * log(1 / 9 * a^2) / (pi * a * nu - pi * a) + 27 / 80 / (pi * a * nu - pi * a)
    f[5, 2, 1] = -3 / 8 * log(25 / 9 * a^2) / (pi * a * nu - pi * a) + 3 / 8 * log(1 / 9 * a^2) / (pi * a * nu - pi * a) + 9 / 8 / (pi * a * nu - pi * a)
    f[5, 3, 1] = 3 / 32 * log(25 / 9 * a^2) / (pi * a * nu - pi * a) - 3 / 32 * log(1 / 9 * a^2) / (pi * a * nu - pi * a) - 9 / 16 / (pi * a * nu - pi * a)
    f[5, 1, 2] = -9 / 16 / (pi * a * nu - pi * a)
    f[5, 2, 2] = 13 / 8 / (pi * a * nu - pi * a)
    f[5, 3, 2] = -9 / 16 / (pi * a * nu - pi * a)
    f[5, 1, 3] = 3 / 32 * log(25 / 9 * a^2) / (pi * a * nu - pi * a) - 3 / 32 * log(1 / 9 * a^2) / (pi * a * nu - pi * a) - 9 / 16 / (pi * a * nu - pi * a)
    f[5, 2, 3] = -3 / 8 * log(25 / 9 * a^2) / (pi * a * nu - pi * a) + 3 / 8 * log(1 / 9 * a^2) / (pi * a * nu - pi * a) + 9 / 8 / (pi * a * nu - pi * a)
    f[5, 3, 3] = 9 / 32 * log(25 / 9 * a^2) / (pi * a * nu - pi * a) - 9 / 32 * log(1 / 9 * a^2) / (pi * a * nu - pi * a) + 27 / 80 / (pi * a * nu - pi * a)

    f[6, 1, 1] = 9 / 32 * log(25 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) - 9 / 32 * log(1 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) + 621 / 100 / (pi * a^2 * nu - pi * a^2)
    f[6, 2, 1] = -9 / 16 * log(25 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) + 9 / 16 * log(1 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) - 27 / 5 / (pi * a^2 * nu - pi * a^2)
    f[6, 3, 1] = 9 / 32 * log(25 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) - 9 / 32 * log(1 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) + 27 / 20 / (pi * a^2 * nu - pi * a^2)
    f[6, 1, 2] = 3 / 4 / (pi * a^2 * nu - pi * a^2)
    f[6, 2, 2] = 0
    f[6, 3, 2] = -3 / 4 / (pi * a^2 * nu - pi * a^2)
    f[6, 1, 3] = -9 / 32 * log(25 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) + 9 / 32 * log(1 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) - 27 / 20 / (pi * a^2 * nu - pi * a^2)
    f[6, 2, 3] = 9 / 16 * log(25 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) - 9 / 16 * log(1 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) + 27 / 5 / (pi * a^2 * nu - pi * a^2)
    f[6, 3, 3] = -9 / 32 * log(25 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) + 9 / 32 * log(1 / 9 * a^2) / (pi * a^2 * nu - pi * a^2) - 621 / 100 / (pi * a^2 * nu - pi * a^2)

    f[7, 1, 1] = -9 / 16 / (a^2 * nu - a^2)
    f[7, 2, 1] = 9 / 8 / (a^2 * nu - a^2)
    f[7, 3, 1] = -9 / 16 / (a^2 * nu - a^2)
    f[7, 1, 2] = -9 / 16 / (a^2 * nu - a^2)
    f[7, 2, 2] = 9 / 8 / (a^2 * nu - a^2)
    f[7, 3, 2] = -9 / 16 / (a^2 * nu - a^2)
    f[7, 1, 3] = -9 / 16 / (a^2 * nu - a^2)
    f[7, 2, 3] = 9 / 8 / (a^2 * nu - a^2)
    f[7, 3, 3] = -9 / 16 / (a^2 * nu - a^2)

    f = permutedims(f, [3, 1, 2]) # This is just a hack until I get around to putting in the correct indices above
    return f
end


export quadkernel_farfield
"""
    quadkernel_farfield(x, y, a, nu)
    
Farfield 3-node quadratic kernels for substitution into CS kernels
f has dimensions [nobs, f=7, shapefunctions=3]

Note: These were calculated symbolically using Maxima.
"""
function quadkernel_farfield(x, y, a, nu)
    f = zeros(length(x), 7, 3)
    for i in 1:length(x)
        # term1 = atan((a - x[i]) / y[i])
        # term2 = atan((a + x[i]) / y[i])
        term1 = pi / 2 * sign(y[i] / (a - x[i])) - atan(y[i] / (a - x[i]))
        term2 = pi / 2 * sign(y[i] / (a + x[i])) - atan(y[i] / (a + x[i]))

        f[i, 1, 1] = -1 / 64 * ( 6 * y[i]^3 * (term2 + term1) - 6 * a^3 * log(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2) - 8 * a^3 - 12 * a^2 * x[i] + 12 * a * x[i]^2 - 12 * a * y[i]^2 + 6 * ( (2 * a * x[i] - 3 * x[i]^2) * term2 + (2 * a * x[i] - 3 * x[i]^2) * term1 ) * y[i] + 3 * (a * x[i]^2 - x[i]^3 - (a - 3 * x[i]) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 3 * (a * x[i]^2 - x[i]^3 - (a - 3 * x[i]) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) ) / (pi * a^2 * nu - pi * a^2)
        f[i, 1, 2] = 1 / 32 * (6 * y[i]^3 * (term2 + term1) + a^3 * log(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2) + a^3 * log(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2) - 8 * a^3 + 12 * a * x[i]^2 - 12 * a * y[i]^2 + 2 * ((4 * a^2 - 9 * x[i]^2) * term2 + (4 * a^2 - 9 * x[i]^2) * term1) * y[i] + (4 * a^2 * x[i] - 3 * x[i]^3 + 9 * x[i] * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (4 * a^2 * x[i] - 3 * x[i]^3 + 9 * x[i] * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^2 * nu - pi * a^2)
        f[i, 1, 3] = -1 / 64 * (6 * y[i]^3 * (term2 + term1) - 6 * a^3 * log(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2) - 8 * a^3 + 12 * a^2 * x[i] + 12 * a * x[i]^2 - 12 * a * y[i]^2 - 6 * ((2 * a * x[i] + 3 * x[i]^2) * term2 + (2 * a * x[i] + 3 * x[i]^2) * term1) * y[i] - 3 * (a * x[i]^2 + x[i]^3 - (a + 3 * x[i]) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + 3 * (a * x[i]^2 + x[i]^3 - (a + 3 * x[i]) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^2 * nu - pi * a^2)

        f[i, 2, 1] = -3 / 32 * (3 * y[i]^2 * (term2 + term1) - (a - 3 * x[i]) * y[i] * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (a - 3 * x[i]) * y[i] * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) - 6 * a * y[i] + (2 * a * x[i] - 3 * x[i]^2) * term2 + (2 * a * x[i] - 3 * x[i]^2) * term1) / (pi * a^2 * nu - pi * a^2)
        f[i, 2, 2] = 1 / 16 * (9 * y[i]^2 * (term2 + term1) + 9 * x[i] * y[i] * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * x[i] * y[i] * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) - 18 * a * y[i] + (4 * a^2 - 9 * x[i]^2) * term2 + (4 * a^2 - 9 * x[i]^2) * term1) / (pi * a^2 * nu - pi * a^2)
        f[i, 2, 3] = -3 / 32 * (3 * y[i]^2 * (term2 + term1) + (a + 3 * x[i]) * y[i] * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (a + 3 * x[i]) * y[i] * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2)) - 6 * a * y[i] - (2 * a * x[i] + 3 * x[i]^2) * term2 - (2 * a * x[i] + 3 * x[i]^2) * term1) / (pi * a^2 * nu - pi * a^2)

        f[i, 3, 1] = 3 / 64 * (8 * a^2 - 12 * a * x[i] - 4 * ((a - 3 * x[i]) * term2 + (a - 3 * x[i]) * term1) * y[i] - (2 * a * x[i] - 3 * x[i]^2 + 3 * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (2 * a * x[i] - 3 * x[i]^2 + 3 * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^2 * nu - pi * a^2)
        f[i, 3, 2] = 1 / 32 * (36 * a * x[i] - 36 * (x[i] * term2 + x[i] * term1) * y[i] + (4 * a^2 - 9 * x[i]^2 + 9 * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (4 * a^2 - 9 * x[i]^2 + 9 * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^2 * nu - pi * a^2)
        f[i, 3, 3] = -3 / 64 * (8 * a^2 + 12 * a * x[i] - 4 * ((a + 3 * x[i]) * term2 + (a + 3 * x[i]) * term1) * y[i] - (2 * a * x[i] + 3 * x[i]^2 - 3 * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (2 * a * x[i] + 3 * x[i]^2 - 3 * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^2 * nu - pi * a^2)

        f[i, 4, 1] = 3 / 32 * (4 * a^2 * y[i]^3 - 2 * ((a - 3 * x[i]) * term2 + (a - 3 * x[i]) * term1) * y[i]^4 - 4 * ((a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * term2 + (a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * term1) * y[i]^2 + 4 * (a^4 - 3 * a^3 * x[i] + a^2 * x[i]^2) * y[i] - 2 * (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5) * term2 - 2 * (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5) * term1 - 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^6 * nu - pi * a^6 + (pi * a^2 * nu - pi * a^2) * x[i]^4 + (pi * a^2 * nu - pi * a^2) * y[i]^4 - 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2 + 2 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^2)
        f[i, 4, 2] = 1 / 16 * (20 * a^3 * x[i] * y[i] - 18 * (x[i] * term2 + x[i] * term1) * y[i]^4 - 36 * ((a^2 * x[i] + x[i]^3) * term2 + (a^2 * x[i] + x[i]^3) * term1) * y[i]^2 - 18 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5) * term2 - 18 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5) * term1 + 9 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^6 * nu - pi * a^6 + (pi * a^2 * nu - pi * a^2) * x[i]^4 + (pi * a^2 * nu - pi * a^2) * y[i]^4 - 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2 + 2 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^2)
        f[i, 4, 3] = -3 / 32 * (4 * a^2 * y[i]^3 - 2 * ((a + 3 * x[i]) * term2 + (a + 3 * x[i]) * term1) * y[i]^4 - 4 * ((a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * term2 + (a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * term1) * y[i]^2 + 4 * (a^4 + 3 * a^3 * x[i] + a^2 * x[i]^2) * y[i] - 2 * (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5) * term2 - 2 * (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5) * term1 + 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 3 * (y[i]^5 + 2 * (a^2 + x[i]^2) * y[i]^3 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * y[i]) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^6 * nu - pi * a^6 + (pi * a^2 * nu - pi * a^2) * x[i]^4 + (pi * a^2 * nu - pi * a^2) * y[i]^4 - 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2 + 2 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^2)

        f[i, 5, 1] = 3 / 32 * (6 * y[i]^5 * (term2 + term1) - 6 * a^5 - 4 * a^4 * x[i] + 18 * a^3 * x[i]^2 + 4 * a^2 * x[i]^3 - 12 * a * x[i]^4 - 12 * a * y[i]^4 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^3 - 2 * (9 * a^3 - 2 * a^2 * x[i] + 12 * a * x[i]^2) * y[i]^2 + 6 * ((a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term2 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term1) * y[i] - (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5 + (a - 3 * x[i]) * y[i]^4 + 2 * (a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + (a^5 - 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 + 6 * a^2 * x[i]^3 + a * x[i]^4 - 3 * x[i]^5 + (a - 3 * x[i]) * y[i]^4 + 2 * (a^3 - 3 * a^2 * x[i] + a * x[i]^2 - 3 * x[i]^3) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^6 * nu - pi * a^6 + (pi * a^2 * nu - pi * a^2) * x[i]^4 + (pi * a^2 * nu - pi * a^2) * y[i]^4 - 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2 + 2 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^2)
        f[i, 5, 2] = -1 / 16 * (18 * y[i]^5 * (term2 + term1) - 26 * a^5 + 62 * a^3 * x[i]^2 - 36 * a * x[i]^4 - 36 * a * y[i]^4 + 36 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^3 - 2 * (31 * a^3 + 36 * a * x[i]^2) * y[i]^2 + 18 * ((a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term2 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term1) * y[i] + 9 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5 + x[i] * y[i]^4 + 2 * (a^2 * x[i] + x[i]^3) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * (a^4 * x[i] - 2 * a^2 * x[i]^3 + x[i]^5 + x[i] * y[i]^4 + 2 * (a^2 * x[i] + x[i]^3) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^6 * nu - pi * a^6 + (pi * a^2 * nu - pi * a^2) * x[i]^4 + (pi * a^2 * nu - pi * a^2) * y[i]^4 - 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2 + 2 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^2)
        f[i, 5, 3] = 3 / 32 * (6 * y[i]^5 * (term2 + term1) - 6 * a^5 + 4 * a^4 * x[i] + 18 * a^3 * x[i]^2 - 4 * a^2 * x[i]^3 - 12 * a * x[i]^4 - 12 * a * y[i]^4 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^3 - 2 * (9 * a^3 + 2 * a^2 * x[i] + 12 * a * x[i]^2) * y[i]^2 + 6 * ((a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term2 + (a^4 - 2 * a^2 * x[i]^2 + x[i]^4) * term1) * y[i] + (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5 + (a + 3 * x[i]) * y[i]^4 + 2 * (a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - (a^5 + 3 * a^4 * x[i] - 2 * a^3 * x[i]^2 - 6 * a^2 * x[i]^3 + a * x[i]^4 + 3 * x[i]^5 + (a + 3 * x[i]) * y[i]^4 + 2 * (a^3 + 3 * a^2 * x[i] + a * x[i]^2 + 3 * x[i]^3) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^6 * nu - pi * a^6 + (pi * a^2 * nu - pi * a^2) * x[i]^4 + (pi * a^2 * nu - pi * a^2) * y[i]^4 - 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2 + 2 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^2)

        f[i, 6, 1] = 3 / 32 * (8 * a^8 - 24 * a^7 * x[i] - 16 * a^6 * x[i]^2 + 60 * a^5 * x[i]^3 + 8 * a^4 * x[i]^4 - 48 * a^3 * x[i]^5 + 12 * a * x[i]^7 + 12 * a * x[i] * y[i]^6 + 4 * (2 * a^4 + 12 * a^3 * x[i] + 9 * a * x[i]^3) * y[i]^4 + 4 * (4 * a^6 + 3 * a^5 * x[i] - 12 * a^4 * x[i]^2 + 9 * a * x[i]^5) * y[i]^2 - 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^10 * nu - pi * a^10 + (pi * a^2 * nu - pi * a^2) * x[i]^8 + (pi * a^2 * nu - pi * a^2) * y[i]^8 - 4 * (pi * a^4 * nu - pi * a^4) * x[i]^6 + 4 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^6 + 6 * (pi * a^6 * nu - pi * a^6) * x[i]^4 + 2 * (3 * pi * a^6 * nu - 3 * pi * a^6 + 3 * (pi * a^2 * nu - pi * a^2) * x[i]^4 + 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2) * y[i]^4 - 4 * (pi * a^8 * nu - pi * a^8) * x[i]^2 + 4 * (pi * a^8 * nu - pi * a^8 + (pi * a^2 * nu - pi * a^2) * x[i]^6 - (pi * a^4 * nu - pi * a^4) * x[i]^4 - (pi * a^6 * nu - pi * a^6) * x[i]^2) * y[i]^2)
        f[i, 6, 2] = 1 / 16 * (56 * a^7 * x[i] - 148 * a^5 * x[i]^3 + 128 * a^3 * x[i]^5 - 36 * a * x[i]^7 - 36 * a * x[i] * y[i]^6 - 12 * (8 * a^3 * x[i] + 9 * a * x[i]^3) * y[i]^4 - 4 * (a^5 * x[i] - 8 * a^3 * x[i]^3 + 27 * a * x[i]^5) * y[i]^2 + 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^10 * nu - pi * a^10 + (pi * a^2 * nu - pi * a^2) * x[i]^8 + (pi * a^2 * nu - pi * a^2) * y[i]^8 - 4 * (pi * a^4 * nu - pi * a^4) * x[i]^6 + 4 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^6 + 6 * (pi * a^6 * nu - pi * a^6) * x[i]^4 + 2 * (3 * pi * a^6 * nu - 3 * pi * a^6 + 3 * (pi * a^2 * nu - pi * a^2) * x[i]^4 + 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2) * y[i]^4 - 4 * (pi * a^8 * nu - pi * a^8) * x[i]^2 + 4 * (pi * a^8 * nu - pi * a^8 + (pi * a^2 * nu - pi * a^2) * x[i]^6 - (pi * a^4 * nu - pi * a^4) * x[i]^4 - (pi * a^6 * nu - pi * a^6) * x[i]^2) * y[i]^2)
        f[i, 6, 3] = -3 / 32 * (8 * a^8 + 24 * a^7 * x[i] - 16 * a^6 * x[i]^2 - 60 * a^5 * x[i]^3 + 8 * a^4 * x[i]^4 + 48 * a^3 * x[i]^5 - 12 * a * x[i]^7 - 12 * a * x[i] * y[i]^6 + 4 * (2 * a^4 - 12 * a^3 * x[i] - 9 * a * x[i]^3) * y[i]^4 + 4 * (4 * a^6 - 3 * a^5 * x[i] - 12 * a^4 * x[i]^2 - 9 * a * x[i]^5) * y[i]^2 + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 + 2 * a * x[i] + x[i]^2 + y[i]^2)) - 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8 + y[i]^8 + 4 * (a^2 + x[i]^2) * y[i]^6 + 2 * (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * y[i]^4 + 4 * (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * y[i]^2) * log(abs(a^2 - 2 * a * x[i] + x[i]^2 + y[i]^2))) / (pi * a^10 * nu - pi * a^10 + (pi * a^2 * nu - pi * a^2) * x[i]^8 + (pi * a^2 * nu - pi * a^2) * y[i]^8 - 4 * (pi * a^4 * nu - pi * a^4) * x[i]^6 + 4 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^6 + 6 * (pi * a^6 * nu - pi * a^6) * x[i]^4 + 2 * (3 * pi * a^6 * nu - 3 * pi * a^6 + 3 * (pi * a^2 * nu - pi * a^2) * x[i]^4 + 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2) * y[i]^4 - 4 * (pi * a^8 * nu - pi * a^8) * x[i]^2 + 4 * (pi * a^8 * nu - pi * a^8 + (pi * a^2 * nu - pi * a^2) * x[i]^6 - (pi * a^4 * nu - pi * a^4) * x[i]^4 - (pi * a^6 * nu - pi * a^6) * x[i]^2) * y[i]^2)

        f[i, 7, 1] = -3 / 16 * (3 * y[i]^8 * (term2 + term1) - 6 * a * y[i]^7 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^6 - 6 * (4 * a^3 + 3 * a * x[i]^2) * y[i]^5 + 6 * ((3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term2 + (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term1) * y[i]^4 - 2 * (15 * a^5 - 8 * a^4 * x[i] + 9 * a * x[i]^4) * y[i]^3 + 12 * ((a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term2 + (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term1) * y[i]^2 - 2 * (6 * a^7 - 8 * a^6 * x[i] + 3 * a^5 * x[i]^2 + 8 * a^4 * x[i]^3 - 12 * a^3 * x[i]^4 + 3 * a * x[i]^6) * y[i] + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term2 + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term1) / (pi * a^10 * nu - pi * a^10 + (pi * a^2 * nu - pi * a^2) * x[i]^8 + (pi * a^2 * nu - pi * a^2) * y[i]^8 - 4 * (pi * a^4 * nu - pi * a^4) * x[i]^6 + 4 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^6 + 6 * (pi * a^6 * nu - pi * a^6) * x[i]^4 + 2 * (3 * pi * a^6 * nu - 3 * pi * a^6 + 3 * (pi * a^2 * nu - pi * a^2) * x[i]^4 + 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2) * y[i]^4 - 4 * (pi * a^8 * nu - pi * a^8) * x[i]^2 + 4 * (pi * a^8 * nu - pi * a^8 + (pi * a^2 * nu - pi * a^2) * x[i]^6 - (pi * a^4 * nu - pi * a^4) * x[i]^4 - (pi * a^6 * nu - pi * a^6) * x[i]^2) * y[i]^2)
        f[i, 7, 2] = 1 / 8 * (9 * y[i]^8 * (term2 + term1) - 18 * a * y[i]^7 + 36 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^6 - 2 * (32 * a^3 + 27 * a * x[i]^2) * y[i]^5 + 18 * ((3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term2 + (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term1) * y[i]^4 - 2 * (37 * a^5 + 8 * a^3 * x[i]^2 + 27 * a * x[i]^4) * y[i]^3 + 36 * ((a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term2 + (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term1) * y[i]^2 - 2 * (14 * a^7 + a^5 * x[i]^2 - 24 * a^3 * x[i]^4 + 9 * a * x[i]^6) * y[i] + 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term2 + 9 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term1) / (pi * a^10 * nu - pi * a^10 + (pi * a^2 * nu - pi * a^2) * x[i]^8 + (pi * a^2 * nu - pi * a^2) * y[i]^8 - 4 * (pi * a^4 * nu - pi * a^4) * x[i]^6 + 4 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^6 + 6 * (pi * a^6 * nu - pi * a^6) * x[i]^4 + 2 * (3 * pi * a^6 * nu - 3 * pi * a^6 + 3 * (pi * a^2 * nu - pi * a^2) * x[i]^4 + 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2) * y[i]^4 - 4 * (pi * a^8 * nu - pi * a^8) * x[i]^2 + 4 * (pi * a^8 * nu - pi * a^8 + (pi * a^2 * nu - pi * a^2) * x[i]^6 - (pi * a^4 * nu - pi * a^4) * x[i]^4 - (pi * a^6 * nu - pi * a^6) * x[i]^2) * y[i]^2)
        f[i, 7, 3] = -3 / 16 * (3 * y[i]^8 * (term2 + term1) - 6 * a * y[i]^7 + 12 * ((a^2 + x[i]^2) * term2 + (a^2 + x[i]^2) * term1) * y[i]^6 - 6 * (4 * a^3 + 3 * a * x[i]^2) * y[i]^5 + 6 * ((3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term2 + (3 * a^4 + 2 * a^2 * x[i]^2 + 3 * x[i]^4) * term1) * y[i]^4 - 2 * (15 * a^5 + 8 * a^4 * x[i] + 9 * a * x[i]^4) * y[i]^3 + 12 * ((a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term2 + (a^6 - a^4 * x[i]^2 - a^2 * x[i]^4 + x[i]^6) * term1) * y[i]^2 - 2 * (6 * a^7 + 8 * a^6 * x[i] + 3 * a^5 * x[i]^2 - 8 * a^4 * x[i]^3 - 12 * a^3 * x[i]^4 + 3 * a * x[i]^6) * y[i] + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term2 + 3 * (a^8 - 4 * a^6 * x[i]^2 + 6 * a^4 * x[i]^4 - 4 * a^2 * x[i]^6 + x[i]^8) * term1) / (pi * a^10 * nu - pi * a^10 + (pi * a^2 * nu - pi * a^2) * x[i]^8 + (pi * a^2 * nu - pi * a^2) * y[i]^8 - 4 * (pi * a^4 * nu - pi * a^4) * x[i]^6 + 4 * (pi * a^4 * nu - pi * a^4 + (pi * a^2 * nu - pi * a^2) * x[i]^2) * y[i]^6 + 6 * (pi * a^6 * nu - pi * a^6) * x[i]^4 + 2 * (3 * pi * a^6 * nu - 3 * pi * a^6 + 3 * (pi * a^2 * nu - pi * a^2) * x[i]^4 + 2 * (pi * a^4 * nu - pi * a^4) * x[i]^2) * y[i]^4 - 4 * (pi * a^8 * nu - pi * a^8) * x[i]^2 + 4 * (pi * a^8 * nu - pi * a^8 + (pi * a^2 * nu - pi * a^2) * x[i]^6 - (pi * a^4 * nu - pi * a^4) * x[i]^4 - (pi * a^6 * nu - pi * a^6) * x[i]^2) * y[i]^2)
    end
    return f
end


export trac2dispstress
"""
    trac2dispstress(xcomp, ycomp, f, y, mu, nu)
    
Calculate interior displacements and stresses from boundary tractions
"""
function trac2dispstress(xcomp, ycomp, f, y, mu, nu)
    disp = zeros(length(y), 2)
    stress = zeros(length(y), 3)
    # _xcomp, _ycomp = -xcomp, -ycomp # For Okada consistency...causes problems...don't do this.
    _xcomp, _ycomp = xcomp, ycomp # No sign flip.  This provides a different standard than Okada but makes more sense for BEM
    # TODO: Check sigma_yy and sigma_xy for forcing in the x-direction
    for i in 1:length(y)
        disp[i, 1] = _xcomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f[i, 1] + y[i] * f[i, 2]) + _ycomp / (2.0 * mu) * (-y[i] * f[i, 3])
        disp[i, 2] = _xcomp / (2.0 * mu) * (-y[i] * f[i, 3]) + _ycomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f[i, 1] - y[i] * f[i, 2])
        stress[i, 1] = _xcomp * ((3.0 - 2.0 * nu) * f[i, 3] + y[i] * f[i, 4]) + _ycomp * (2.0 * nu * f[i, 2] + y[i] * -f[i, 5])
        # stress[i, 2] = _xcomp * (-1.0 * (1.0 - 2.0 * nu) * f[i, 3] + y[i] * f[i, 4]) + _ycomp * (2.0 * (1.0 - nu) * f[i, 2] - y[i] * -f[i, 5])
        # stress[i, 3] = _xcomp * (2.0 * (1.0 - nu) * f[i, 2] + y[i] * f[i, 5]) + _ycomp * ((1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4])
        stress[i, 2] = _xcomp * (-1.0 * (1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4]) + _ycomp * (2.0 * (1.0 - nu) * f[i, 2] - y[i] * -f[i, 5])
        stress[i, 3] = _xcomp * (2.0 * (1.0 - nu) * f[i, 2] + y[i] * -f[i, 5]) + _ycomp * ((1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4])
    end
    return disp, stress
end


export slip2dispstress
"""
    slip2dispstress(xcomp, ycomp, f, y, mu, nu)
    
Calculate interior displacements and stresses from boundary slip/displacement
"""
function slip2dispstress(xcomp, ycomp, f, y, mu, nu)
    disp, stress = zeros(length(y), 2), zeros(length(y), 3)
    _xcomp, _ycomp = -xcomp, -ycomp # For Okada consistency
    # _xcomp, _ycomp = xcomp, ycomp # TODO: Neglecting the sign flip leads to errors with ...super strange
    for i in 1:length(y)
        # CS equations (5.2.3), page 81
        disp[i, 1] = _xcomp * (2.0 * (1.0 - nu) * f[i, 2] - y[i] * f[i, 5]) +
                     _ycomp * (-(1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4])
        disp[i, 2] = _xcomp * ((1.0 - 2.0 * nu) * f[i, 3] - y[i] * f[i, 4]) +
                     _ycomp * (2.0 * (1 - nu) * f[i, 2] - y[i] * -f[i, 5])

        # CS equations (5.2.4), page 81
        stress[i, 1] = 2.0 * _xcomp * mu * (2.0 * f[i, 4] + y[i] * f[i, 6]) +
                       2.0 * _ycomp * mu * (-f[i, 5] + y[i] * f[i, 7])
        stress[i, 2] = 2.0 * _xcomp * mu * (-y[i] * f[i, 6]) +
                       2.0 * _ycomp * mu * (-f[i, 5] - y[i] * f[i, 7])
        stress[i, 3] = 2.0 * _xcomp * mu * (-f[i, 5] + y[i] * f[i, 7]) +
                       2.0 * _ycomp * mu * (-y[i] * f[i, 6])
    end
    return disp, stress
end


export stress2trac
"""
    stress2trac(stress, nvec)
    
Calculate traction from stresses and normal vector.
"""
function stress2trac(stress, nvec)
    return [stress[1] stress[3] ; stress[3] stress[2]] * nvec
end


export standardize_elements!
"""
    standardize_elements!(els)

In place updating els structure to determin derived parameters
(e.g., midpoint, length, etc.)

Note: Should be called after declaring elements and before any 
BEM calculations (boundary or interior)
"""
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


export rotdispstress
"""
    rotdispstress(disp, stress, rotmat)

Rotate displacement vectors and stress tensors.

Note: Consider moving to StaticArrays
"""
function rotdispstress(disp, stress, rotmat)
    # TODO: If this is slow hand expand the matrix vector mulltiplies
    # TODO: inplace for speed.  Some benchmarks suggest 50x speedup!
    _disp = zeros(size(disp))
    _stress = zeros(size(stress))
    for i in 1:size(stress)[1]
        _disp[i, 1], _disp[i, 2] = rotmat * [disp[i, 1] ; disp[i, 2]]
        stresstensor = [stress[i, 1] stress[i, 3] ; stress[i, 3] stress[i, 2]]
        _stress[i, 1], _stress[i, 3], _, _stress[i, 2] = rotmat * stresstensor * transpose(rotmat)
    end
    return _disp, _stress
end


export plotelements
"""
    plotelements(els)

Plot all non-NaN parts of elements as straing lines.
"""
function plotelements(els)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0)
    end
    return nothing
end


"""
    stylesubplots(xlim, ylim)

Common styling for subplots generated by plotfields.
"""
function stylesubplots(xlim, ylim)
    gca().set_aspect("equal")
    gca().set_xlim([xlim[1], xlim[2]])
    gca().set_ylim([ylim[1], ylim[2]])
    gca().set_xticks([xlim[1], xlim[2]])
    gca().set_yticks([ylim[1], ylim[2]])
    return nothing
end


"""
    plotfields_contours(els, xobs, yobs, idx, field, title)

Contouring of a single field for six-panel plotfields.
"""
function plotfields_contours(els, xobs, yobs, idx, field, title)
    ncontours = 40
    xlim = [minimum(xobs) maximum(xobs)]
    ylim = [minimum(yobs) maximum(yobs)]
    subplot(2, 3, idx)
    scale = 1.0
    fieldmax = maximum(@.abs(field))
    contourf(xobs, yobs, reshape(field, size(xobs)), ncontours,
        vmin = -scale * fieldmax, vmax = scale * fieldmax, cmap = rycroftcmap())
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
"""
    plotfields(els, xobs, yobs, disp, stress, suptitle)

Set up six panel plots for displacement and stress fields.
"""
function plotfields(els, xobs, yobs, disp, stress, suptitle)
    figure(figsize = (25, 15))
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
    PyPlot.suptitle(suptitle, fontsize=20)
    show()
    return nothing
end


export rycroftcmap
"""
    rycroftcmap()

Chris Rycroft designed and shared this colormap.
"""
function rycroftcmap()
    cmap = ColorScheme([Colors.RGB(255.0 / 255.0, 20.0 / 255.0, 40.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 47.0 / 255.0, 146.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 138.0 / 255.0, 216.0 / 255.0),
                        Colors.RGB(55.0 / 255.0, 145.0 / 255.0, 230.0 / 255.0),
                        Colors.RGB(150.0 / 255.0, 230.0 / 255.0, 80.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 251.0 / 255.0, 0.0 / 255.0),
                        Colors.RGB(255.0 / 255.0, 255.0 / 255.0, 255.0 / 255.0)],
                        "rycroft", "psychedelic")
    cmap = ColorMap(ColorScheme([get(cmap, i) for i in 0.0:0.001:1.0]).colors)
    return cmap
end


export PUTC
"""
    PUTC(fun2dispstress, els, obsidx, srcidx, mu, nu)

Convenience function to calcualtion partial derivates of displacement, 
and tractions for constant slip/displacement elements.  This 
simply swaps the order of arguments (srcidx, obsidx) -> (obsidx, srcidx)
and does not return the stress partials.
"""
function PUTC(fun2dispstress, els, obsidx, srcidx, mu, nu)
    T, _, H = partialsconstdispstress(fun2dispstress, els, srcidx, obsidx, mu, nu)
    return T, H
end


export PUTQ
"""
    PUTQ(fun2dispstress, els, obsidx, srcidx, mu, nu)

Conveinece function to calcualtion partial derivates of displacement, 
and tractions for constant slip/displacement elements.  This 
simply swaps the order of arguments (srcidx, obsidx) -> (obsidx, srcidx)
and does not return the stress partials.
"""
function PUTQ(fun2dispstress, els, obsidx, srcidx, mu, nu)
    T, _, H = partialsquaddispstress(fun2dispstress, els, srcidx, obsidx, mu, nu)
    return T, H
end

    
export partialsconstdispstress
"""
    partialsconstdispstress(fun2dispstress, els, srcidx, obsidx, mu, nu)

Partial derivates of displacement, stress and tractions for constant
slip/displacement elements.

Note: Rename to parUSTC?
"""
function partialsconstdispstress(fun2dispstress, els, srcidx, obsidx, mu, nu)
    nobs, nsrc = length(obsidx), length(srcidx)
    partialsdisp, partialsstress, partialstrac = zeros(2 * nobs, 2 * nsrc), zeros(3 * nobs, 2 * nsrc), zeros(2 * nobs, 2 * nsrc)
    _partialsdisp, _partialsstress, _partialstrac = zeros(2, 2), zeros(3, 2), zeros(2, 2)
    for isrc in 1:nsrc
        for iobs in 1:nobs
            _partialsdisp[:, 1], _partialsstress[:, 1] = constdispstress(fun2dispstress, els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]], els, srcidx[isrc], 1, 0, mu, nu)
            _partialsdisp[:, 2], _partialsstress[:, 2] = constdispstress(fun2dispstress, els.xcenter[obsidx[iobs]], els.ycenter[obsidx[iobs]], els, srcidx[isrc], 0, 1, mu, nu)
            _partialstrac[:, 1] = stress2trac(_partialsstress[:, 1], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            _partialstrac[:, 2] = stress2trac(_partialsstress[:, 2], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            partialsdisp[2 * (iobs - 1) + 1:2 * (iobs - 1) + 2, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _partialsdisp
            partialsstress[3 * (iobs - 1) + 1:3 * (iobs - 1) + 3, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _partialsstress
            partialstrac[2 * (iobs - 1) + 1:2 * (iobs - 1) + 2, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _partialstrac
        end
    end
    return partialsdisp, partialsstress, partialstrac
end


export quaddispstressinterleaved
"""
    quaddispstressinterleaved(fun2dispstress, xobs, yobs, els, idx, para, perp, mu, nu)

This is small helper function to return interleaved displacements
and stresses rather than the stacked columns.  For the non-coincident case.
"""
function quaddispstressinterleaved(fun2dispstress, xobs, yobs, els, idx, para, perp, mu, nu)
    dispstacked, stressstacked = quaddispstress(fun2dispstress, xobs, yobs, els, idx, para, perp, mu, nu)
    dispinterleaved = zeros(2 * length(xobs))
    dispinterleaved[1:2:end] = dispstacked[:, 1]
    dispinterleaved[2:2:end] = dispstacked[:, 2]
    stressinterleaved = zeros(3 * length(xobs))
    stressinterleaved[1:3:end] = stressstacked[:, 1]
    stressinterleaved[2:3:end] = stressstacked[:, 2]
    stressinterleaved[3:3:end] = stressstacked[:, 3]
    return dispinterleaved, stressinterleaved
end


export quaddispstresscoincidentinterleaved
"""
    quaddispstresscoincidentinterleaved(fun2dispstress, xobs, yobs, els, idx, para, perp, mu, nu)

This is small helper function to return interleaved displacements
and stresses rather than the stacked columns.  For the coincident case.
"""
function quaddispstresscoincidentinterleaved(fun2dispstress, xobs, yobs, els, idx, para, perp, mu, nu)
    dispstacked, stressstacked = quaddispstresscoincident(fun2dispstress, xobs, yobs, els, idx, para, perp, mu, nu)
    dispinterleaved = zeros(2 * length(xobs))
    dispinterleaved[1:2:end] = dispstacked[:, 1]
    dispinterleaved[2:2:end] = dispstacked[:, 2]
    stressinterleaved = zeros(3 * length(xobs))
    stressinterleaved[1:3:end] = stressstacked[:, 1]
    stressinterleaved[2:3:end] = stressstacked[:, 2]
    stressinterleaved[3:3:end] = stressstacked[:, 3]
    return dispinterleaved, stressinterleaved
end


export pointonline
"""
    pointonline(x1, y1, x2, y2, x3)

Given two endpointts of a line and an x value somewhere in between
calculate and return the corresponding y value
"""
function pointonline(x1, y1, x2, y2, x3)
    y3 = ((y2 - y1) * (x3 - x2)) / (x2 - x1) + y2
    return y3
end


export partialsquaddispstress
"""
    partialsquaddispstress(fun2dispstress, els, srcidx, obsidx, mu, nu)

Calculate partials derivaties for quadratic elements.

Each 6x6 displacement submatrix is shaped as:
       sx(ϕ1) sy(ϕ1) sx(ϕ2) sy(ϕ2) sx(ϕ3) sy(ϕ3)
ux(ϕ1)
uy(ϕ1)
ux(ϕ2)
uy(ϕ2)
ux(ϕ3)
uy(ϕ3)

Note: Rename to: parUSTQ?
"""
function partialsquaddispstress(fun2dispstress, els, srcidx, obsidx, mu, nu)
    nobs = length(obsidx)
    nsrc = length(srcidx)
    partialsdisp = zeros(6 * nobs, 6 * nsrc)
    partialsstress = zeros(9 * nobs, 6 * nsrc)
    partialstrac = zeros(6 * nobs, 6 * nsrc)
    _partialsdisp = zeros(6, 6)
    _partialsstress = zeros(9, 6)
    _partialstrac = zeros(6, 6)

    # Far-field interactions for all except main (block) diagonal
    for isrc in 1:nsrc
        for iobs in 1:nobs
            if srcidx[isrc] == obsidx[iobs]
                _partialsdisp[:, 1], _partialsstress[:, 1] = quaddispstresscoincidentinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [1 0 0], [0 0 0], mu, nu)
                _partialsdisp[:, 2], _partialsstress[:, 2] = quaddispstresscoincidentinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [1 0 0], mu, nu)
                _partialsdisp[:, 3], _partialsstress[:, 3] = quaddispstresscoincidentinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 1 0], [0 0 0], mu, nu)
                _partialsdisp[:, 4], _partialsstress[:, 4] = quaddispstresscoincidentinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 1 0], mu, nu)
                _partialsdisp[:, 5], _partialsstress[:, 5] = quaddispstresscoincidentinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 1], [0 0 0], mu, nu)
                _partialsdisp[:, 6], _partialsstress[:, 6] = quaddispstresscoincidentinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 0 1], mu, nu)
            else
                _partialsdisp[:, 1], _partialsstress[:, 1] = quaddispstressinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [1 0 0], [0 0 0], mu, nu)
                _partialsdisp[:, 2], _partialsstress[:, 2] = quaddispstressinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [1 0 0], mu, nu)
                _partialsdisp[:, 3], _partialsstress[:, 3] = quaddispstressinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 1 0], [0 0 0], mu, nu)
                _partialsdisp[:, 4], _partialsstress[:, 4] = quaddispstressinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 1 0], mu, nu)
                _partialsdisp[:, 5], _partialsstress[:, 5] = quaddispstressinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 1], [0 0 0], mu, nu)
                _partialsdisp[:, 6], _partialsstress[:, 6] = quaddispstressinterleaved(fun2dispstress, els.xnodes[obsidx[iobs], :], els.ynodes[obsidx[iobs], :], els, srcidx[isrc], [0 0 0], [0 0 1], mu, nu)
            end
            for i in 1:6
                _partialstrac[1:2, i] = stress2trac(_partialsstress[1:3, i], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
                _partialstrac[3:4, i] = stress2trac(_partialsstress[4:6, i], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
                _partialstrac[5:6, i] = stress2trac(_partialsstress[7:9, i], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            end
            partialsdisp[6 * (iobs - 1) + 1:6 * (iobs - 1) + 6, 6 * (isrc - 1) + 1:6 * (isrc - 1) + 6] = _partialsdisp
            partialsstress[9 * (iobs - 1) + 1:9 * (iobs - 1) + 9, 6 * (isrc - 1) + 1:6 * (isrc - 1) + 6] = _partialsstress
            partialstrac[6 * (iobs - 1) + 1:6 * (iobs - 1) + 6, 6 * (isrc - 1) + 1:6 * (isrc - 1) + 6] = _partialstrac
        end
    end
    return partialsdisp, partialsstress, partialstrac
end


"""
    slip2coef(xnodes, slip, a)

Given a set of 3 displacements/slips at quadratic nodes convert
to a set of three phi values neccesary to comput the three quadratic
functions.
"""
function slip2coef(xnodes, slip, a)
    mat = zeros(3, 3)
    mat[:, 1] = (xnodes ./ a) .* (9 .* (xnodes ./ a) ./ 8 .- 3 ./ 4)
    mat[:, 2] = (1 .- 3 .* (xnodes ./ a) ./ 2) .* (1 .+ 3 .* (xnodes ./ a) ./ 2)
    mat[:, 3] = (xnodes ./ a) .* (9 .* (xnodes ./ a) ./ 8 .+ 3 ./ 4)
    return inv(mat) * slip
end


"""
    coef2slip(xnodes, slip, a)

Given a set of 3 phi values at quadratic nodes convert
to a set of three displacement/slip values neccesary to 
physically interpret the results.
"""
function coef2slip(xnodes, coef, a)
    mat = zeros(3, 3)
    mat[:, 1] = (xnodes ./ a) .* (9 .* (xnodes ./ a) ./ 8 .- 3 ./ 4)
    mat[:, 2] = (1 .- 3 .* (xnodes ./ a) ./ 2) .* (1 .+ 3 * (xnodes ./ a) ./ 2)
    mat[:, 3] = (xnodes ./ a) .* (9 .* (xnodes ./ a) ./ 8 .+ 3 ./ 4)
    return mat * coef
end


export getidx
"""
    getidx(label, els)

Return a vector with the indices of els.name that match label.
"""
function getidx(label, els)
    idx = findall(x->x == label, els.name)
    println("getidx found " * string(length(idx)) * " elements with label \"" * label * "\"")
    return idx
end


export getidxdict
"""
    getidxdict(els)

Return a dictionary of the arrays containing the indices of els
associated with each distinct name in els.  The keys to the 
dictionary are each unique name in els.
"""
function getidxdict(els)
    idxdict = Dict()
    names = unique(els.name)
    for i in 1:length(names)
        if length(names[i]) > 0
            idxdict[names[i]] = getidx(names[i], els)
        end
    end
    return idxdict
end


export initpartials
"""
    initpartials(els)

Create a nested dictionary for storing partial derivatives.
Should probably change this to Hstar, Sstar, and Tstar notation.
"""
function initpartials(els)
    partials = Dict()
    partials["disp"] = Dict()
    partials["stress"] = Dict()
    partials["trac"] = Dict()
    fieldnames = collect(keys(partials))
    elnames = collect(keys(getidxdict(els)))
    for iname in 1:length(fieldnames)
        for i in 1:length(elnames)
            partials[fieldnames[iname]][elnames[i]] = Dict()
            for j in 1:length(elnames)
                partials[fieldnames[iname]][elnames[i]][elnames[j]] = []
            end
        end
    end
    return partials
end


export quadstack
"""
    quadstack(vec)

Utility function for "stacking" a flattened quad vector
Useful for converting results of BEM solve to arguments
that can be used with quaddispstress for forward models
"""
function quadstack(vec)
    stack = transpose(reshape(vec, 3, Int(length(vec) / 3)))
    return stack
end


export plotqdtimeseries
"""
    plotqdtimeseries(sol, stridesize, nels)

Utility function for plotting quasidynaic results. For constant slip elements.
Can run out of memory :(
"""
function plotqdtimeseries(sol, stridesize, nels)
    siay = 365.25 * 24 * 60 * 60
    t = [x / siay for x in sol.t]
    theta = zeros(length(t), nels)
    vx = zeros(length(t), nels)
    vy = zeros(length(t), nels)
    for i in 1:length(t)
        vx[i, :] = sol.u[i][1:stridesize:end]
        vy[i, :] = sol.u[i][2:stridesize:end]
        theta[i, :] = sol.u[i][3:stridesize:end]
    end

    figure(1, figsize = (15, 8))
    subplot(3, 2, 1)
    plot(t, vx, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_x")

    subplot(3, 2, 3)
    plot(t, vy, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_y")

    subplot(3, 2, 5)
    plot(t, theta, "-", linewidth = 0.5)
    yscale("log")
    xlabel("t (years)")
    ylabel(L"\theta")

    subplot(3, 2, 2)
    plot(1:1:length(t), vx, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_x")

    subplot(3, 2, 4)
    plot(1:1:length(t), vy, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_y")

    subplot(3, 2, 6)
    plot(1:1:length(t), theta, "-", linewidth = 0.5)
    yscale("log")
    xlabel("time step #")
    ylabel(L"\theta")

    figure(2, figsize = (15, 5))
    plotme = log10.(vx')
    contourf(plotme, 5, cmap = rycroftcmap())
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}v$ (m/s)")
    contour(plotme, 5, linewidths = 0.5, linestyles = "solid", colors = "k")
    xlabel("time step")
    ylabel("element index")
    show()
end


export plotqdtimeseriesquad
"""
    plotqdtimeseriesquad(sol, stridesize, nels)

Utility function for plotting quasidynaic results. For quadratic slip elements.
Can run out of memory :(
"""
function plotqdtimeseriesquad(sol, stridesize, nels)
    siay = 365.25 * 24 * 60 * 60
    t = [x / siay for x in sol.t]
    theta = zeros(length(t), 3 * nels)
    vx = zeros(length(t), 3 * nels)
    vy = zeros(length(t), 3 * nels)
    for i in 1:length(t)
        vx[i, :] = sol.u[i][1:stridesize:end]
        vy[i, :] = sol.u[i][2:stridesize:end]
        theta[i, :] = sol.u[i][3:stridesize:end]
    end

    figure(figsize = (15, 8))
    subplot(3, 2, 1)
    plot(t, vx, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_x")

    subplot(3, 2, 3)
    plot(t, vy, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_y")

    subplot(3, 2, 5)
    plot(t, theta, "-", linewidth = 0.5)
    yscale("log")
    xlabel("t (years)")
    ylabel(L"\theta")

    subplot(3, 2, 2)
    plot(1:1:length(t), vx, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_x")

    subplot(3, 2, 4)
    plot(1:1:length(t), vy, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_y")

    subplot(3, 2, 6)
    plot(1:1:length(t), theta, "-", linewidth = 0.5)
    yscale("log")
    xlabel("time step #")
    ylabel(L"\theta")

    figure(figsize = (15, 5))
    plotme = log10.(vx')
    contourf(plotme, 5, cmap = rycroftcmap())
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}v$ (m/s)")
    contour(plotme, 5, linewidths = 0.5, linestyles = "solid", colors = "k")
    xlabel("time step")
    ylabel("element index")
    show()
end


export thetaaginglaw
"""
    thetaaginglaw(v, theta, dc)

Aging law in rate and state friction
"""
function thetaaginglaw(v, theta, dc)
    return @. 1 - theta * v / dc
end


export thetasliplaw
"""
    thetasliplaw(v, theta, dc)

Slip evolution law in rate and state friction
"""
function thetasliplaw(v, theta, dc)
    return @. -v * theta / dc * log(v * theta / dc)
end


export gravityparticularfunctions
"""
    gravityparticularfunctions()

From Pape and Bannerjee 1987
"""
function gravityparticularfunctions(x, y, g, rho, lambda, mu)
    U = zeros(length(x), 2)
    S = zeros(length(x), 3)
    U[:, 1] = @. -lambda * rho * g / (4 * mu * (lambda + mu)) * x * y  # Pape and Banerjee (1987) equation (6a)
    U[:, 2] = @. (rho * g) / (8 * mu * (lambda + mu)) * (lambda * x^2 + (lambda + 2 * mu) * y^2) # Pape and Banerjee (1987) equation (6b)
    S[:, 1] .= 0 # Pape and Banerjee (1987) equation (6c)
    S[:, 2] = @. rho * g * y  # Pape and Banerjee (1987) equation (6d)
    S[:, 3] .= 0 # Pape and Banerjee (1987) equation (6e)
    return U, S
end


end
