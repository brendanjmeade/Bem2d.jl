module Bem2d

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


# def disp_stres_const(
#     x,
#     y,
#     a,
#     mu,
#     nu,
#     shape_function,
#     element_type,
#     x_component,
#     y_component,
#     x_center,
#     y_center,
#     rotation_matrix,
#     inverse_rotation_matrix,
# ):
#     """ Calculate displacements and stresses for constant and linear slip elements """
#     # Rotate and translate into local coordinate system
#     x = x - x_center
#     y = y - y_center
#     rotated_coords = np.matmul(np.vstack((x, y)).T, rotation_matrix)
#     x = rotated_coords[:, 0]
#     y = rotated_coords[:, 1]
#
#     # Convert to global coordinates here.  Should this be elsewhere?
#     global_components = inverse_rotation_matrix @ np.array([x_component, y_component])
#     x_component = global_components[0]
#     y_component = global_components[1]
#
#     if shape_function == "constant":
#         f = constant_kernel(x, y, a, nu)
#     elif shape_function == "linear":
#         f = linear_kernel(x, y, a, nu)
#
#     if element_type == "traction":
#         displacement, stress = f_traction_to_displacement_stress(
#             x_component, y_component, f, y, mu, nu
#         )
#     elif element_type == "slip":
#         displacement, stress = f_slip_to_displacement_stress(
#             x_component, y_component, f, y, mu, nu
#         )
#
#     displacement, stress = rotate_displacement_stress(
#         displacement, stress, inverse_rotation_matrix
#     )
#     return displacement, stress
#
#
# def constant_kernel(x, y, a, nu):
#     """ From Starfield and Crouch, pages 49 and 82 """
#     f = np.zeros((7, x.size))
#
#     f[0, :] = (
#         -1
#         / (4 * np.pi * (1 - nu))
#         * (
#             y * (np.arctan2(y, (x - a)) - np.arctan2(y, (x + a)))
#             - (x - a) * np.log(np.sqrt((x - a) ** 2 + y ** 2))
#             + (x + a) * np.log(np.sqrt((x + a) ** 2 + y ** 2))
#         )
#     )
#
#     f[1, :] = (
#         -1
#         / (4 * np.pi * (1 - nu))
#         * ((np.arctan2(y, (x - a)) - np.arctan2(y, (x + a))))
#     )
#
#     f[2, :] = (
#         1
#         / (4 * np.pi * (1 - nu))
#         * (
#             np.log(np.sqrt((x - a) ** 2 + y ** 2))
#             - np.log(np.sqrt((x + a) ** 2 + y ** 2))
#         )
#     )
#
#     f[3, :] = (
#         1
#         / (4 * np.pi * (1 - nu))
#         * (y / ((x - a) ** 2 + y ** 2) - y / ((x + a) ** 2 + y ** 2))
#     )
#
#     f[4, :] = (
#         1
#         / (4 * np.pi * (1 - nu))
#         * ((x - a) / ((x - a) ** 2 + y ** 2) - (x + a) / ((x + a) ** 2 + y ** 2))
#     )
#
#     f[5, :] = (
#         1
#         / (4 * np.pi * (1 - nu))
#         * (
#             ((x - a) ** 2 - y ** 2) / ((x - a) ** 2 + y ** 2) ** 2
#             - ((x + a) ** 2 - y ** 2) / ((x + a) ** 2 + y ** 2) ** 2
#         )
#     )
#
#     f[6, :] = (
#         2
#         * y
#         / (4 * np.pi * (1 - nu))
#         * (
#             (x - a) / ((x - a) ** 2 + y ** 2) ** 2
#             - (x + a) / ((x + a) ** 2 + y ** 2) ** 2
#         )
#     )
#     return f

# def f_traction_to_displacement_stress(x_component, y_component, f, y, mu, nu):
#     """ This is the generalization from Starfield and Crouch """
#     displacement = np.zeros((2, y.size))
#     stress = np.zeros((3, y.size))
#
#     # The sign change here is to:
#     # 1 - Ensure consistenty with Okada convention
#     # 2 - For a horizontal/flat fault make the upper half move in the +x direction
#     x_component = -1 * x_component
#     y_component = -1 * y_component
#
#     displacement[0, :] = x_component / (2 * mu) * (
#         (3 - 4 * nu) * f[0, :] + y * f[1, :]
#     ) + y_component / (2 * mu) * (-y * f[2, :])
#
#     displacement[1, :] = x_component / (2 * mu) * (-y * f[2, :]) + y_component / (
#         2 * mu
#     ) * ((3 - 4 * nu) * f[0, :] - y * f[1, :])
#
#     stress[0, :] = x_component * (
#         (3 - 2 * nu) * f[2, :] + y * f[3, :]
#     ) + y_component * (2 * nu * f[1, :] + y * f[4, :])
#
#     stress[1, :] = x_component * (
#         -1 * (1 - 2 * nu) * f[2, :] + y * f[3, :]
#     ) + y_component * (2 * (1 - nu) * f[1, :] - y * f[4, :])
#
#     stress[2, :] = x_component * (
#         2 * (1 - nu) * f[1, :] + y * f[4, :]
#     ) + y_component * ((1 - 2 * nu) * f[2, :] - y * f[3, :])
#
#     return displacement, stress
# end

# Generalization from Starfield and Crouch
export slip2dispstress
function slip2dispstress(xcomp, ycomp, f, y, mu, nu)
    disp = zeros(length(y), 2)
    stress = zeros(length(y), 3)
    _xcomp = -xcomp # For Okada consistency
    _ycomp = -ycomp # For Okada consistency
    for i in 1:length(y)
        disp[i, 1] = _xcomp * (2.0 * (1.0 - nu) * f[i, 2] - y * f[i, 5]) + _ycomp * (-1.0 * (1.0 - 2.0 * nu) * f[i, 3] - y * f[i, 4])
        disp[i, 2] = _xcomp * ((1.0 - 2.0 * nu) * f[i, 3] - y * f[i, 4]) + _ycomp * (2.0 * (1 - nu) * f[i, 2] - y * -f[i, 5])
        stress[i, 1] = 2.0 * _xcomp * mu * (2.0 * f[i, 4] + y * f[i, 6]) + 2.0 * _ycomp * mu * (-f[i, 5] + y * f[i, 7])
        stress[i, 2] = 2.0 * _xcomp * mu * (-y * f[i, 6]) + 2.0 * _ycomp * mu * (-f[i, 5] - y * f[i, 7])
        stress[i, 3] = 2.0 * _xcomp * mu * (-f[i, 5] + y * f[i, 7]) + 2.0 * _ycomp * mu * (-y * f[i, 6])
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
function rotdispstress(disp, stress, rotmat_inv)
    # If this is slow hand expand the matrix vector multiplies
    # inplace for speed.  Some benchmarks suggest 50x speedup!
    for i in 0:size(stress)[1]
        dispglobal = rotmat_inv * [disp[i, 1] ; disp[i, 2]]
        _disp[i, 1] = dispglobal[1]
        _disp[i, 2] = dispglobal[2]
        stress_tensor = [stress[i, 1] stress[i, 3] ; stress[i, 3] stress[i, 2]]
        stressglobal = rotmat_inv' * stress_tensor * rotmat_inv
        _stress[i, 1] = stressglobal[1, 1]
        _stress[i, 2] = stressglobal[2, 2]
        _stress[i, 3] = stressglobal[1, 2]
    end
    return _disp, _stress
end

end
