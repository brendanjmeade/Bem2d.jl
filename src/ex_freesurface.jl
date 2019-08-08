import numpy as np
import matplotlib.pyplot as plt
import bem2d
from importlib import reload
from okada_wrapper import dc3dwrapper

bem2d = reload(bem2d)

plt.close("all")

# Material properties and observation grid
SURFACE_SIGN_FLIP = True
mu = 30e9
nu = 0.25
n_pts = 102
width = 5
x = np.linspace(-width, width, n_pts)
y = np.linspace(-width, width, n_pts)
x, y = np.meshgrid(x, y)
x = x.flatten()
y = y.flatten()

# Define elements
elements_surface = []
elements_fault = []
element = {}

# Traction free surface
x1, y1, x2, y2 = bem2d.discretized_line(-5, 0, 5, 0, 20)
for i in range(0, x1.size):
    element["x1"] = x1[i]
    element["y1"] = y1[i]
    element["x2"] = x2[i]
    element["y2"] = y2[i]
    element["name"] = "surface"
    elements_surface.append(element.copy())
elements_surface = bem2d.standardize_elements(elements_surface)

# Constant slip fault
x1, y1, x2, y2 = bem2d.discretized_line(-1, -1, 0, 0, 1)
for i in range(0, x1.size):
    element["x1"] = x1[i]
    element["y1"] = y1[i]
    element["x2"] = x2[i]
    element["y2"] = y2[i]
    element["name"] = "fault"
    elements_fault.append(element.copy())
elements_fault = bem2d.standardize_elements(elements_fault)

d1, s1, t1 = bem2d.constant_partials_all(elements_surface, elements_fault, mu, nu)
d2, s2, t2 = bem2d.constant_partials_all(elements_surface, elements_surface, mu, nu)
d1_quadratic, s1_quadratic, t1_quadratic = bem2d.quadratic_partials_all(
    elements_surface, elements_fault, mu, nu
)
d2_quadratic, s2_quadratic, t2_quadratic = bem2d.quadratic_partials_all(
    elements_surface, elements_surface, mu, nu
)

# Constant case: Predict surface displacements from unit strike slip forcing
x_center = np.array([_["x_center"] for _ in elements_surface])
fault_slip = np.zeros(2 * len(elements_fault))
fault_slip[0::2] = np.sqrt(2) / 2
fault_slip[1::2] = np.sqrt(2) / 2
disp_full_space = d1 @ fault_slip
disp_free_surface = np.linalg.inv(t2) @ (t1 @ fault_slip)

# Quadratic case: Predict surface displacements from unit strike slip forcing
x_center_quadratic = np.array(
    [_["x_integration_points"] for _ in elements_surface]
).flatten()
fault_slip_quadratic = np.zeros(6 * len(elements_fault))
fault_slip_quadratic[0::2] = np.sqrt(2) / 2
fault_slip_quadratic[1::2] = np.sqrt(2) / 2
disp_full_space_quadratic = d1_quadratic @ fault_slip_quadratic
disp_free_surface_quadratic = np.linalg.inv(t2_quadratic) @ (
    t1_quadratic @ fault_slip_quadratic
)

# Okada solution for 45 degree dipping fault
x_okada = np.linspace(-5, 5, 1000)
disp_okada_x = np.zeros(x_okada.shape)
disp_okada_y = np.zeros(x_okada.shape)

for i in range(0, x_okada.size):
    # Fault dipping at 45 degrees
    _, u, _ = dc3dwrapper(
        2.0 / 3.0,
        [0, x_okada[i] + 0.5, 0],
        0.5,
        45,  # 135
        [-1000, 1000],
        [-np.sqrt(2) / 2, np.sqrt(2) / 2],
        [0.0, 1.0, 0.0],
    )
    disp_okada_x[i] = u[1]
    disp_okada_y[i] = u[2]

plt.figure(figsize=(6, 8))
plt.subplot(2, 1, 1)
plt.plot(x_okada, disp_okada_x, "-k", linewidth=0.5, label="Okada")
plt.plot(
    x_center,
    disp_free_surface[0::2],
    "bo",
    markerfacecolor="None",
    linewidth=0.5,
    label="constant BEM",
)
plt.plot(
    x_center_quadratic,
    disp_free_surface_quadratic[0::2],
    "r.",
    linewidth=0.5,
    label="quadratic BEM",
)
plt.xlim([-5, 5])
plt.ylim([-1, 1])
plt.xticks(np.arange(-5, 6))
plt.yticks(np.linspace(-1, 1, 9))
plt.xlabel(r"$x$")
plt.ylabel("displacement")
plt.title(r"$u_x$")
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(x_okada, disp_okada_y, "-k", linewidth=0.5, label="Okada")
plt.plot(
    x_center,
    disp_free_surface[1::2],
    "bo",
    markerfacecolor="None",
    linewidth=0.1,
    label="constant BEM",
)
plt.plot(
    x_center_quadratic,
    disp_free_surface_quadratic[1::2],
    "r.",
    linewidth=0.5,
    label="quadratic BEM",
)

plt.xlim([-5, 5])
plt.ylim([-1, 1])
plt.xticks(np.arange(-5, 6))
plt.yticks(np.linspace(-1, 1, 9))
plt.xlabel(r"$x$")
plt.ylabel("displacement")
plt.title(r"$u_y$")
plt.legend()
plt.tight_layout()


# Okada internal dispalcements
displacement_okada = np.zeros((2, x.size))
stress_okada = np.zeros((3, x.size))
disp_okada_x = np.zeros(x.shape)
disp_okada_y = np.zeros(x.shape)
stress_okada_xx = np.zeros(x.shape)
stress_okada_yy = np.zeros(x.shape)
stress_okada_xy = np.zeros(x.shape)

for i in range(0, x.size):
    _, u, s = dc3dwrapper(
        2.0 / 3.0,
        [0, x[i] + 0.5, y[i]],
        0.5,
        45,  # 135
        [-1000, 1000],
        [-np.sqrt(2) / 2, np.sqrt(2) / 2],
        [0.0, 1.0, 0.0],
    )
    disp_okada_x[i] = u[1]
    disp_okada_y[i] = u[2]
    dgt_xx = s[1, 1]
    dgt_yy = s[2, 2]
    dgt_xy = s[1, 2]
    dgt_yx = s[2, 1]
    e_xx = dgt_xx
    e_yy = dgt_yy
    e_xy = 0.5 * (dgt_yx + dgt_xy)
    s_xx = mu * (e_xx + e_yy) + 2 * mu * e_xx
    s_yy = mu * (e_xx + e_yy) + 2 * mu * e_yy
    s_xy = 2 * mu * e_xy
    stress_okada_xx[i] = s_xx
    stress_okada_yy[i] = s_yy
    stress_okada_xy[i] = s_xy

displacement_okada[0, :] = disp_okada_x
displacement_okada[1, :] = disp_okada_y
stress_okada[0, :] = stress_okada_xx
stress_okada[1, :] = stress_okada_yy
stress_okada[2, :] = stress_okada_xy

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_okada,
    stress_okada,
    "Okada",
)

#
# Internal evaluation for constant BEM
#
fault_slip_ss = fault_slip[0::2]
fault_slip_ts = fault_slip[1::2]
displacement_full_space = np.zeros((2, x.size))
stress_full_space = np.zeros((3, x.size))
for i, element in enumerate(elements_fault):
    displacement, stress = bem2d.displacements_stresses_constant_linear(
        x,
        y,
        element["half_length"],
        mu,
        nu,
        "constant",
        "slip",
        fault_slip_ss[i].copy(),
        fault_slip_ts[i].copy(),
        element["x_center"],
        element["y_center"],
        element["rotation_matrix"],
        element["inverse_rotation_matrix"],
    )
    displacement_full_space += displacement
    stress_full_space += stress

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_full_space,
    stress_full_space,
    "constant BEM: fault full space",
)

# Free surface
fault_slip_x = disp_free_surface[0::2]
fault_slip_y = disp_free_surface[1::2]
displacement_free_surface = np.zeros((2, x.size))
stress_free_surface = np.zeros((3, x.size))
for i, element in enumerate(elements_surface):
    displacement, stress = bem2d.displacements_stresses_constant_linear(
        x,
        y,
        element["half_length"],
        mu,
        nu,
        "constant",
        "slip",
        fault_slip_x[i].copy(),
        fault_slip_y[i].copy(),
        element["x_center"],
        element["y_center"],
        element["rotation_matrix"],
        element["inverse_rotation_matrix"],
    )
    displacement_free_surface += displacement
    stress_free_surface += stress

if SURFACE_SIGN_FLIP:
    displacement_free_surface *= -1
    stress_free_surface *= -1

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_free_surface,
    stress_free_surface,
    "constant BEM: free surface",
)

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_free_surface + displacement_full_space,
    stress_free_surface + stress_full_space,
    "constant BEM: fault + free surface",
)

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_free_surface + displacement_full_space - displacement_okada,
    stress_free_surface + stress_full_space - stress_okada,
    "constant BEM - Okada",
)


#
# Internal evaluation for Quadratic BEM
#
fault_slip_ss_quadratic = fault_slip_quadratic[0::2]
fault_slip_ts_quadratic = fault_slip_quadratic[1::2]
displacement_quadratic_elements = np.zeros((2, x.size))
stress_quadratic_elements = np.zeros((3, x.size))
for i, element in enumerate(elements_fault):
    displacement, stress = bem2d.displacements_stresses_quadratic_NEW(
        x,
        y,
        element["half_length"],
        mu,
        nu,
        "slip",
        fault_slip_ss_quadratic[i * 3 : (i + 1) * 3].copy(),
        fault_slip_ts_quadratic[i * 3 : (i + 1) * 3].copy(),
        element["x_center"],
        element["y_center"],
        element["rotation_matrix"],
        element["inverse_rotation_matrix"],
    )
    displacement_quadratic_elements += displacement
    stress_quadratic_elements += stress

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_quadratic_elements,
    stress_quadratic_elements,
    "quadratic BEM: fault full space",
)

# Free surface
surface_slip_x_quadratic = disp_free_surface_quadratic[0::2]
surface_slip_y_quadratic = disp_free_surface_quadratic[1::2]
displacement_free_surface_quadratic = np.zeros((2, x.size))
stress_free_surface_quadratic = np.zeros((3, x.size))
for i, element in enumerate(elements_surface):
    displacement, stress = bem2d.displacements_stresses_quadratic_NEW(
        x,
        y,
        element["half_length"],
        mu,
        nu,
        "slip",
        surface_slip_x_quadratic[i * 3 : (i + 1) * 3].copy(),
        surface_slip_y_quadratic[i * 3 : (i + 1) * 3].copy(),
        element["x_center"],
        element["y_center"],
        element["rotation_matrix"],
        element["inverse_rotation_matrix"],
    )
    displacement_free_surface_quadratic += displacement
    stress_free_surface_quadratic += stress

if SURFACE_SIGN_FLIP:
    displacement_free_surface_quadratic *= -1
    stress_free_surface_quadratic *= -1

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_free_surface_quadratic,
    stress_free_surface_quadratic,
    "quadratic BEM: free surface",
)

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_free_surface_quadratic + displacement_quadratic_elements,
    stress_free_surface_quadratic + stress_quadratic_elements,
    "quadratic BEM: fault + free surface",
)

bem2d.plot_fields(
    elements_surface + elements_fault,
    x.reshape(n_pts, n_pts),
    y.reshape(n_pts, n_pts),
    displacement_free_surface_quadratic
    + displacement_quadratic_elements
    - displacement_okada,
    stress_free_surface_quadratic + stress_quadratic_elements - stress_okada,
    "quadratic BEM - Okada",
)
