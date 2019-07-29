module Bem2d

export Elements
mutable struct Elements
    x1::Array{Float64, 1}
    y1::Array{Float64, 1}
    x2::Array{Float64, 1}
    y2::Array{Float64, 1}
    name::Array{String, 1}
    uxglobal::Array{Float64, 1}
    uyglobal::Array{Float64, 1}
    angle::Array{Float64, 1}
    length::Array{Float64, 1}
    halflength::Array{Float64, 1}
    xcenter::Array{Float64, 1}
    ycenter::Array{Float64, 1}
    rotationmatrix::Array{Float64, 1}
    rotationmatrixinverse::Array{Float64, 1}
    xnormal::Array{Float64, 1}
    ynormal::Array{Float64, 1}
    xnodes::Array{Float64, 1}
    ynodes::Array{Float64, 1}
    Elements() = new([], [], [], [], [], [],
                     [], [], [], [], [], [],
                     [], [], [], [], [], [])

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


# def displacements_stresses_constant_linear(
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
#
#
# def f_slip_to_displacement_stress(x_component, y_component, f, y, mu, nu):
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
#     displacement[0, :] = x_component * (
#         2 * (1 - nu) * f[1, :] - y * f[4, :]
#     ) + y_component * (-1 * (1 - 2 * nu) * f[2, :] - y * f[3, :])
#
#     displacement[1, :] = x_component * (
#         (1 - 2 * nu) * f[2, :] - y * f[3, :]
#     ) + y_component * (
#         2 * (1 - nu) * f[1, :] - y * -f[4, :]
#     )  # Note the negative sign in front f[4, :] because f[4, :] = f,xx = -f,yy
#
#     stress[0, :] = 2 * x_component * mu * (
#         2 * f[3, :] + y * f[5, :]
#     ) + 2 * y_component * mu * (-f[4, :] + y * f[6, :])
#
#     stress[1, :] = 2 * x_component * mu * (-y * f[5, :]) + 2 * y_component * mu * (
#         -f[4, :] - y * f[6, :]
#     )
#
#     stress[2, :] = 2 * x_component * mu * (
#         -f[4, :] + y * f[6, :]
#     ) + 2 * y_component * mu * (-y * f[5, :])
#
#     return displacement, stress


# def stress_to_traction(stress, normal_vector):
#     """ Compute tractions from stress tensor and normal vector """
#     stress_tensor = np.zeros((2, 2))
#     stress_tensor[0, 0] = stress[0]
#     stress_tensor[0, 1] = stress[2]
#     stress_tensor[1, 0] = stress[2]
#     stress_tensor[1, 1] = stress[1]
#     traction = stress_tensor @ normal_vector
#     return traction

export standardize_elements!
function standardize_elements!(elements)
    for i in 1:length(elements.x1)
        dx = elements.x2[i] - elements.x1[i]
        dy = elements.y2[i] - elements.y1[i]
        mag = sqrt(dx^2 + dy^2)
        push!(elements.angle, atan(elements.y2[i] - elements.y1[i], elements.x2[i] - elements.x1[i]))
        push!(elements.length, sqrt((elements.x2[i] - elements.x1[i])^2 + (elements.y2[i] - elements.y1[i])^2))
        push!(elements.halflength, 0.5 * elements.length[i])
        push!(elements.xcenter, 0.5 * (elements.x2[i] + elements.x1[i]))
        push!(elements.ycenter, 0.5 * (elements.y2[i] + elements.y1[i]))
        # element["rotation_matrix"] = np.array(
        #     [
        #         [np.cos(element["angle"]), -np.sin(element["angle"])],
        #         [np.sin(element["angle"]), np.cos(element["angle"])],
        #     ]
        # )
        # element["inverse_rotation_matrix"] = np.array(
        #     [
        #         [np.cos(-element["angle"]), -np.sin(-element["angle"])],
        #         [np.sin(-element["angle"]), np.cos(-element["angle"])],
        #     ]
        # )
        push!(elements.xnormal, dy / mag)
        push!(elements.ynormal, -dx / mag)
        elements.xnodes = vcat(elements.xnodes, [elements.xcenter[i] - (2 / 3 * dx / 2), elements.xcenter[i], elements.xcenter[i] + (2 / 3 * dx / 2)])
        # element["ynodes"] = np.array(
        #     [
        #         element["y_center"] - (2 / 3 * dy / 2),
        #         element["y_center"],
        #         element["y_center"] + (2 / 3 * dy / 2),
        #     ]
        # )
        #
        # # If a local boundary condition is giving convert to global
        # # TODO: This is just for convenience there should be flags for real BCs
        # if "ux_local" in element:
        #     u_local = np.array([element["ux_local"], element["uy_local"]])
        #     u_global = element["rotation_matrix"] @ u_local
        #     element["ux_global_constant"] = u_global[0]
        #     element["uy_global_constant"] = u_global[1]
        #     element["ux_global_quadratic"] = np.repeat(u_global[0], 3)
        #     element["uy_global_quadratic"] = np.repeat(u_global[1], 3)
    end
    return nothing
end

# def rotate_displacement_stress(displacement, stress, inverse_rotation_matrix):
#     """ Rotate displacements stresses from local to global reference frame """
#     displacement = np.matmul(displacement.T, inverse_rotation_matrix).T
#     for i in range(0, stress.shape[1]):
#         stress_tensor = np.array(
#             [[stress[0, i], stress[2, i]], [stress[2, i], stress[1, i]]]
#         )
#         stress_tensor_global = (
#             inverse_rotation_matrix.T @ stress_tensor @ inverse_rotation_matrix
#         )
#         stress[0, i] = stress_tensor_global[0, 0]
#         stress[1, i] = stress_tensor_global[1, 1]
#         stress[2, i] = stress_tensor_global[0, 1]
#     return displacement, stress



end
