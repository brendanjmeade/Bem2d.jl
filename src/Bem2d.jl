module Bem2d

# http://julia-programming-language.2336112.n4.nabble.com/Meshgrid-function-td37003.html
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
    lastidx::Int64
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
                0) # lastidx
end

export updatelastidx!
function updatelastidx!(elements)
    elements.lastidx = findall(isnan, elements.x1)[1] - 1
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
    # Rotate and translate into local coordinate system with *global* slip components
    _x = zeros(length(x))
    _y = zeros(length(y))
    for i in 1:length(x)
        _x[i], _y[i] = rotmat * [x[i] - xcenter ; y[i] - ycenter]
    end
    _xcomp, _ycomp = rotmatinv * [xcomp ; ycomp]
    f = constantkernel(_x, _y, a, nu)
    disp, stress = slip2dispstress(_xcomp, _ycomp, f, _y, mu, nu)
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
export traction2dispstress
function traction2dispstress(xcomp, ycomp, f, y, mu, nu)
    disp = zeros(length(y), 2)
    stress = zeros(length(y), 3)
    _xcomp = -xcomp # For Okada consistency
    _ycomp = -ycomp # For Okada consistency
    disp[:, 1] = _xcomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f[:, 1] + y * f[:, 2]) + _ycomp / (2.0 * mu) * (-y * f[:, 3])
    disp[:, 2] = _xcomp / (2.0 * mu) * (-y * f[:, 3]) + _ycomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f[:, 1] - y * f[:, 2])
    stress[:, 1] = _xcomp * ((3.0 - 2.0 * nu) * f[:, 3] + y * f[:, 4]) + _ycomp * (2.0 * nu * f[:, 2] + y * f[:, 5])
    stress[:, 2] = _xcomp * (-1.0 * (1.0 - 2.0 * nu) * f[:, 3] + y * f[:, 4]) + _ycomp * (2.0 * (1.0 - nu) * f[:, 2] - y * f[:, 5])
    stress[:, 3] = _xcomp * (2.0 * (1.0 - nu) * f[:, 2] + y * f[:, 5]) + _ycomp * ((1.0 - 2.0 * nu) * f[:, 3] - y * f[:, 4])
    return disp, stress
end

# Generalization from Starfield and Crouch
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

export stress2traction
function stress2traction(stress, nvec)
    traction = [stress[1] stress[2] ; stress[2] stress[3]] * nvec
    return traction
end

export standardize_elements!
function standardize_elements!(elements)
    updatelastidx!(elements)
    for i in 1:elements.lastidx
        dx = elements.x2[i] - elements.x1[i]
        dy = elements.y2[i] - elements.y1[i]
        magnitude = sqrt(dx^2 + dy^2)
        elements.angle[i] = atan(dy, dx)
        elements.length[i] = magnitude
        elements.halflength[i] = 0.5 * elements.length[i]
        elements.xcenter[i] = 0.5 * (elements.x2[i] + elements.x1[i])
        elements.ycenter[i] = 0.5 * (elements.y2[i] + elements.y1[i])
        elements.rotmat[i, :, :] = [cos(elements.angle[1]) -sin(elements.angle[1]) ; sin(elements.angle[1]) cos(elements.angle[1])]
        elements.rotmatinv[i, :, :] = [cos(-elements.angle[1]) -sin(-elements.angle[1]) ; sin(-elements.angle[1]) cos(-elements.angle[1])]
        elements.xnormal[i] = dy / magnitude
        elements.ynormal[i] = -dx / magnitude
        elements.xnodes[i, :] = [elements.xcenter[i] - (2 / 3 * dx / 2), elements.xcenter[i], elements.xcenter[i] + (2 / 3 * dx / 2)]
        elements.ynodes[i, :] = [elements.ycenter[i] - (2 / 3 * dy / 2), elements.ycenter[i], elements.ycenter[i] + (2 / 3 * dy / 2)]
    end
    return nothing
end

export rotdispstress
function rotdispstress(disp, stress, rotmatinv)
    # TODO: If this is slow hand expand the matrix vector multiplies
    # inplace for speed.  Some benchmarks suggest 50x speedup!
    _disp = zeros(size(disp))
    _stress = zeros(size(stress))
    for i in 1:size(stress)[1]
        _disp[i, 1], _disp[i, 2] = rotmatinv * [disp[i, 1] ; disp[i, 2]]
        stresstensor = [stress[i, 1] stress[i, 3] ; stress[i, 3] stress[i, 2]]
        _stress[i, 1], _stress[i, 3], _stress[i, 3] = rotmatinv' * stresstensor * rotmatinv
    end
    return _disp, _stress
end

end
