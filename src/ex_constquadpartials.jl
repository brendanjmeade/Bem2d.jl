import time
import numpy as np
import pandas as pd
from importlib import reload
import bem2d
import matplotlib.pyplot as plt

bem2d = reload(bem2d)

# Material and geometric constants
plt.close("all")
SLIP_TYPE = "tensile_slip"
print(SLIP_TYPE)
mu = 3e10
nu = 0.25
n_elements = 10
elements = []
element = {}
L = 10000
x1, y1, x2, y2 = bem2d.discretized_line(-L, 0, L, 0, n_elements)

for i in range(0, x1.size):
    element["x1"] = x1[i]
    element["y1"] = y1[i]
    element["x2"] = x2[i]
    element["y2"] = y2[i]
    elements.append(element.copy())
elements = bem2d.standardize_elements(elements)


# Convert to dictionary of lists
elements_DL = pd.DataFrame(elements).to_dict("list")
# Convert back to list of dictionaries
elements_LD = pd.DataFrame(elements_DL).to_dict("records")
print(elements == elements_LD)

partials_displacement_constant, partials_stress_constant, partials_traction_constant = bem2d.constant_partials_all(
    elements, elements, mu, nu
)

partials_displacement_quadratic, partials_stress_quadratic, partials_traction_quadratic = bem2d.quadratic_partials_all(
    elements, elements, mu, nu
)

# Create fault elements
N_ELEMENTS = 10
N_NODES = 3 * N_ELEMENTS
ELEMENTS_FAULT = []
ELEMENT = {}
L = 10000
x1, y1, x2, y2 = bem2d.discretized_line(-L, 0, L, 0, N_ELEMENTS)

for i in range(0, x1.size):
    ELEMENT["x1"] = x1[i]
    ELEMENT["y1"] = y1[i]
    ELEMENT["x2"] = x2[i]
    ELEMENT["y2"] = y2[i]
    ELEMENT["a"] = 0.015  # frictional a parameter ()
    ELEMENT["b"] = 0.020  # frictional b parameter ()
    ELEMENT["sigma_n"] = 50e6  # normal stress (Pa)
    ELEMENTS_FAULT.append(ELEMENT.copy())
ELEMENTS_FAULT = bem2d.standardize_elements(ELEMENTS_FAULT)

PARTIALS_DISPLACEMENT_CONSTANT, PARTIALS_STRESS_CONSTANT, PARTIALS_TRACTION_CONSTANT = bem2d.constant_partials_all(
    ELEMENTS_FAULT, ELEMENTS_FAULT, mu, nu
)

PARTIALS_DISPLACEMENT_QUADRATIC, PARTIALS_STRESS_QUADRATIC, PARTIALS_TRACTION_QUADRATIC = bem2d.quadratic_partials_all(
    ELEMENTS_FAULT, ELEMENTS_FAULT, mu, nu
)

# Evaluation points and slip
x_eval = np.array([_["x_integration_points"] for _ in elements]).flatten()
y_eval = np.array([_["y_integration_points"] for _ in elements]).flatten()
slip_quadratic = np.zeros(6 * n_elements)
slip_constant = np.zeros(2 * n_elements)

# Strike-slip
if SLIP_TYPE == "strike_slip":
    slip_quadratic[0::2] = 1  # constant x-slip global
    slip_constant[0::2] = 1  # constant x-slip global
elif SLIP_TYPE == "tensile_slip":
    slip_quadratic[1::2] = 1  # constant x-slip global
    slip_constant[1::2] = 1  # constant x-slip global

# Predict displacements, stresses, and tractions
displacement_quadratic = partials_displacement_quadratic @ slip_quadratic
stress_quadratic = partials_stress_quadratic @ slip_quadratic
traction_quadratic = partials_traction_quadratic @ slip_quadratic
displacement_constant = partials_displacement_constant @ slip_constant
stress_constant = partials_stress_constant @ slip_constant
traction_constant = partials_traction_constant @ slip_constant

# Plot geometry of elements
plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1)
for element in elements:
    plt.plot(
        [element["x1"], element["x2"]],
        [element["y1"], element["y2"]],
        "-k",
        color="r",
        linewidth=0.5,
    )
    plt.plot(
        [element["x1"], element["x2"]],
        [element["y1"], element["y2"]],
        "r.",
        markersize=1,
        linewidth=0.5,
    )
for i, element in enumerate(elements):
    plt.text(
        element["x_center"],
        element["y_center"],
        str(i),
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=8,
    )
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("element geometry")

# Plot predicted displacements
plt.subplot(2, 2, 2)
plt.plot(
    x_eval[1::3],
    displacement_quadratic[2::6],
    "+r",
    markeredgewidth=0.5,
    label="u_x quadratic",
)
plt.plot(
    x_eval[1::3],
    displacement_quadratic[3::6],
    "+b",
    markeredgewidth=0.5,
    label="u_y quadratic",
)
plt.plot(
    x_eval[1::3],
    displacement_constant[0::2],
    "or",
    markerfacecolor="none",
    markeredgewidth=0.5,
    label="u_x constant",
)
plt.plot(
    x_eval[1::3],
    displacement_constant[1::2],
    "ob",
    markerfacecolor="none",
    markeredgewidth=0.5,
    label="u_y constant",
)
plt.legend(loc="upper right")
plt.xlabel("x (m)")
plt.ylabel("displacements (m)")
plt.title("displacements")

plt.subplot(2, 2, 3)
plt.plot(x_eval, slip_quadratic[0::2], "+k", label="quadratic", markeredgewidth=0.5)
plt.plot(
    x_eval[1::3],
    slip_constant[0::2],
    "ok",
    markerfacecolor="none",
    label="constant",
    markeredgewidth=0.5,
)
plt.legend()
plt.xlabel("x (m)")
plt.ylabel("input slip")
plt.title("element slip")

plt.subplot(2, 2, 4)
plt.plot(
    x_eval, stress_quadratic[0::3], "+r", label="s_xx quadratic", markeredgewidth=0.5
)
plt.plot(
    x_eval, stress_quadratic[1::3], "+b", label="s_yy quadratic", markeredgewidth=0.5
)
plt.plot(
    x_eval, stress_quadratic[2::3], "+k", label="s_xy quadratic", markeredgewidth=0.5
)

plt.plot(
    x_eval[1::3],
    stress_constant[0::3],
    "or",
    markerfacecolor="none",
    markeredgewidth=0.5,
    label="s_xx constant",
)
plt.plot(
    x_eval[1::3],
    stress_constant[1::3],
    "ob",
    markerfacecolor="none",
    markeredgewidth=0.5,
    label="s_yy constant",
)
plt.plot(
    x_eval[1::3],
    stress_constant[2::3],
    "ok",
    markerfacecolor="none",
    markeredgewidth=0.5,
    label="s_xy constant",
)
plt.legend()
plt.xlabel("x (m)")
plt.ylabel("stresses (Pa)")
plt.title("stresses")
plt.show(block=False)
