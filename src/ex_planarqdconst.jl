from importlib import reload
import copy
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.integrate

import cppimport.import_hook
import bem2d
import bem2d.newton_rate_state

bem2d = reload(bem2d)
plt.close("all")

# TODO: Try inclined fault
# TODO: Figure out how to enter boundary condition better?

# Constants and model parameters
OUTDIR = "/Users/meade/Desktop/output/"
NEWTON_TOL = 1e-12
MAXITER = 5000
ODE_ATOL = 1e-6
ODE_RTOL = 1e-6
N_NODES_PER_ELEMENT = 3
SPY = 365.25 * 24 * 60 * 60  # Seconds per year
TIME_INTERVAL = SPY * np.array([0.0, 1000.0])
PARAMETERS = {}
PARAMETERS["mu"] = 3e10
PARAMETERS["nu"] = 0.25
PARAMETERS["density"] = 2700  # rock density (kg/m^3)
PARAMETERS["cs"] = np.sqrt(
    PARAMETERS["mu"] / PARAMETERS["density"]
)  # Shear wave speed (m/s)
PARAMETERS["eta"] = PARAMETERS["mu"] / (
    2 * PARAMETERS["cs"]
)  # Radiation damping coefficient (kg / (m^2 * s))
PARAMETERS["d_c"] = 0.05  # state evolution length scale (m)
PARAMETERS["f_0"] = 0.6  # baseline coefficient of friction
PARAMETERS["block_velocity_x"] = 1e-9
PARAMETERS["block_velocity_y"] = 0
PARAMETERS["v_0"] = 1e-6  # reference velocity

# Create fault elements
N_ELEMENTS = 50
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

# Extract parameters for later use as arrays
ELEMENTS_FAULT_ARRAYS = {}
ELEMENTS_FAULT_ARRAYS["a"] = np.array(
    [np.tile(_["a"], 3) for _ in ELEMENTS_FAULT]
).flatten()
ELEMENTS_FAULT_ARRAYS["b"] = np.array(
    [np.tile(_["b"], 3) for _ in ELEMENTS_FAULT]
).flatten()
ELEMENTS_FAULT_ARRAYS["additional_normal_stress"] = np.array(
    [np.tile(_["sigma_n"], 3) for _ in ELEMENTS_FAULT]
).flatten()
ELEMENTS_FAULT_ARRAYS["element_normals"] = np.empty((N_ELEMENTS, 2))
ELEMENTS_FAULT_ARRAYS["element_normals"][:, 0] = np.array(
    [_["x_normal"] for _ in ELEMENTS_FAULT]
).flatten()
ELEMENTS_FAULT_ARRAYS["element_normals"][:, 1] = np.array(
    [_["y_normal"] for _ in ELEMENTS_FAULT]
).flatten()

# Calculate slip to traction partials on the fault
START_TIME = time.time()
_, SLIP_TO_STRESS, SLIP_TO_TRACTION = bem2d.quadratic_partials_all(
    ELEMENTS_FAULT, ELEMENTS_FAULT, PARAMETERS["mu"], PARAMETERS["nu"]
)
PARTIALS_TIME = time.time() - START_TIME


def calc_state(state, velocity):
    """ State evolution law : aging law """
    dstate_dt = (
        ELEMENTS_FAULT_ARRAYS["b"]
        * PARAMETERS["v_0"]
        / PARAMETERS["d_c"]
        * (
            np.exp((PARAMETERS["f_0"] - state) / ELEMENTS_FAULT_ARRAYS["b"])
            - (velocity / PARAMETERS["v_0"])
        )
    )
    return dstate_dt


def steady_state(velocities):  # Should I eliminate this since its only called once?
    """ Steady-state state for initial condition """
    state = scipy.optimize.fsolve(calc_state, 0.5 * np.ones(N_NODES), args=(velocities))
    return state


def calc_derivatives(t, x_and_state):
    """ Derivatives to feed to ODE integrator """
    state = x_and_state[2::3]
    x = np.empty(2 * N_NODES)
    x[0::2] = x_and_state[0::3]
    x[1::2] = x_and_state[1::3]

    # Current shear stress on fault (slip->traction)
    tractions = SLIP_TO_TRACTION @ x

    # Solve for the current velocity...This is the algebraic part
    current_velocity = np.empty(2 * N_NODES)
    bem2d.newton_rate_state.rate_state_solver(
        ELEMENTS_FAULT_ARRAYS["element_normals"],
        tractions,
        state,
        current_velocity,  # Modified in place and returns velocity solution
        ELEMENTS_FAULT_ARRAYS["a"],
        PARAMETERS["eta"],
        PARAMETERS["v_0"],
        0.0,
        ELEMENTS_FAULT_ARRAYS["additional_normal_stress"],
        NEWTON_TOL,
        MAXITER,
        N_NODES_PER_ELEMENT,
    )

    dx_dt = -current_velocity  # Is the negative sign for slip deficit convention?
    dx_dt[0::2] += PARAMETERS["block_velocity_x"]
    dx_dt[1::2] += PARAMETERS["block_velocity_y"]
    velocity_magnitude = np.linalg.norm(current_velocity.reshape((-1, 2)), axis=1)
    dstate_dt = calc_state(state, velocity_magnitude)
    derivatives = np.empty(3 * N_NODES)
    derivatives[0::3] = dx_dt[0::2]
    derivatives[1::3] = dx_dt[1::2]
    derivatives[2::3] = dstate_dt
    return derivatives


# Set initial conditions and time integrate
INITIAL_VELOCITY_X = 1e-3 * PARAMETERS["block_velocity_x"] * np.ones(N_NODES)
INITIAL_VELOCITY_Y = 1e-3 * PARAMETERS["block_velocity_y"] * np.ones(N_NODES)
INITIAL_VELOCITY_MAGNITUDE = np.sqrt(INITIAL_VELOCITY_X ** 2 + INITIAL_VELOCITY_Y ** 2)
INITIAL_CONDITIONS = np.empty(3 * N_NODES)
INITIAL_CONDITIONS[0::3] = INITIAL_VELOCITY_X
INITIAL_CONDITIONS[1::3] = INITIAL_VELOCITY_Y
INITIAL_CONDITIONS[2::3] = steady_state(INITIAL_VELOCITY_MAGNITUDE)

SOLVER = scipy.integrate.RK23(
    calc_derivatives,
    TIME_INTERVAL.min(),
    INITIAL_CONDITIONS,
    TIME_INTERVAL.max(),
    rtol=ODE_RTOL,
    atol=ODE_ATOL,
    vectorized=False,
)

START_TIME = time.time()
SOLUTION = {}
SOLUTION["t"] = [SOLVER.t]
SOLUTION["y"] = [SOLVER.y.copy()]
SOLUTION["stress"] = [np.zeros(3 * N_NODES)]

while SOLVER.t < TIME_INTERVAL.max():
    SOLVER.step()
    N_STEP = len(SOLUTION["t"])
    print(
        f"t = {SOLVER.t / SPY:05.6f}"
        + " of "
        + f"{TIME_INTERVAL.max() / SPY:05.6f}"
        + f" ({100 * SOLVER.t / TIME_INTERVAL.max():.3f}"
        + "%)"
        + f" {N_STEP} steps"
    )
    SOLUTION["t"].append(SOLVER.t)
    SOLUTION["y"].append(SOLVER.y.copy())

    # Calculate on fault stresses too
    x = np.empty(2 * N_NODES)
    x[0::2] = SOLVER.y[0::3]
    x[1::2] = SOLVER.y[1::3]
    stresses = SLIP_TO_STRESS @ x
    SOLUTION["stress"].append(stresses.copy())

SOLUTION["t"] = np.array(SOLUTION["t"])
SOLUTION["y"] = np.array(SOLUTION["y"])
# SOLUTION['u_x'] = SOLUTION['y'][:,0::2]
SOLUTION["stress"] = np.array(SOLUTION["stress"])
SOLVER_TIME = time.time() - START_TIME
print("Partials time: " + f"{PARTIALS_TIME:.2f}" + " (seconds)")
print("Solver time: " + f"{SOLVER_TIME:.2f}" + " (seconds)")


def plot_time_series():
    """ Plot time integrated time series for each node """
    plt.figure(figsize=(12, 9))
    plt.subplot(3, 2, 1)
    plt.plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 0::3], linewidth=0.5)
    plt.ylabel("$u_x$ (m)")

    plt.subplot(3, 2, 2)
    plt.plot(SOLUTION["y"][:, 0::3], linewidth=0.5)
    plt.ylabel("$u_x$ (m)")

    plt.subplot(3, 2, 3)
    plt.plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 1::3], linewidth=0.5)
    plt.ylabel("$u_y$ (m)")

    plt.subplot(3, 2, 4)
    plt.plot(SOLUTION["y"][:, 1::3], linewidth=0.5)
    plt.ylabel("$u_y$ (m)")

    plt.subplot(3, 2, 5)
    plt.plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 2::3], linewidth=0.5)
    plt.xlabel("time (years)")
    plt.ylabel("state")

    plt.subplot(3, 2, 6)
    plt.plot(SOLUTION["y"][:, 2::3], linewidth=0.5)
    plt.xlabel("steps")
    plt.ylabel("state")
    plt.suptitle("Displacement and state")
    plt.show(block=False)


plot_time_series()


def plot_time_series_velocity():
    """ Plot time integrated time series for each node """
    plt.figure(figsize=(12, 9))
    t_diff = SOLUTION["t"][0:-1]  # WHAAA???
    dt = np.diff(SOLUTION["t"])
    solution_diff = np.diff(SOLUTION["y"], axis=0) / dt[:, None]
    solution_diff = np.log10(np.abs(solution_diff))

    plt.subplot(3, 2, 1)
    plt.plot(t_diff / SPY, solution_diff[:, 0::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}v_x$ (m/s)")

    plt.subplot(3, 2, 2)
    plt.plot(solution_diff[:, 0::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}v_x$ (m/s)")

    plt.subplot(3, 2, 3)
    plt.plot(t_diff / SPY, solution_diff[:, 1::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}v_y$ (m/s)")

    plt.subplot(3, 2, 4)
    plt.plot(solution_diff[:, 1::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}v_y$ (m/s)")

    plt.subplot(3, 2, 5)
    plt.plot(t_diff / SPY, solution_diff[:, 2::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}\dot{\theta}$ (1/s)")
    plt.xlabel("$t$ (years)")

    plt.subplot(3, 2, 6)
    plt.plot(solution_diff[:, 2::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}\dot{\theta}$ (1/s)")
    plt.xlabel("step")
    plt.suptitle("Displacement and state rate changes")
    plt.show(block=False)


plot_time_series_velocity()


def plot_stress_time_series():
    """ Plot time integrated time series for each node """
    plt.figure(figsize=(12, 9))
    plt.subplot(3, 2, 1)
    plt.plot(SOLUTION["t"] / SPY, SOLUTION["stress"][:, 0::3], linewidth=0.5)
    plt.ylabel(r"$\sigma_{xx}$ (Pa)")

    plt.subplot(3, 2, 2)
    plt.plot(SOLUTION["stress"][:, 0::3], linewidth=0.5)
    plt.ylabel(r"$\sigma_{xx}$ (Pa)")

    plt.subplot(3, 2, 3)
    plt.plot(SOLUTION["t"] / SPY, SOLUTION["stress"][:, 1::3], linewidth=0.5)
    plt.ylabel(r"$\sigma_{yy}$ (Pa)")

    plt.subplot(3, 2, 4)
    plt.plot(SOLUTION["stress"][:, 1::3], linewidth=0.5)
    plt.ylabel(r"$\sigma_{yy}$ (Pa)")

    plt.subplot(3, 2, 5)
    plt.plot(SOLUTION["t"] / SPY, SOLUTION["stress"][:, 2::3], linewidth=0.5)
    plt.xlabel("time (years)")
    plt.ylabel(r"$\sigma_{xy}$ (Pa)")

    plt.subplot(3, 2, 6)
    plt.plot(SOLUTION["stress"][:, 2::3], linewidth=0.5)
    plt.xlabel("steps")
    plt.ylabel(r"$\sigma_{xy}$ (Pa)")
    plt.suptitle("Stresses")
    plt.show(block=False)


plot_stress_time_series()


def plot_stress_time_series_velocity():
    """ Plot time integrated time series for each node """
    t_diff = SOLUTION["t"][0:-1]  # WHAAA???
    dt = np.diff(SOLUTION["t"])
    solution_diff = np.diff(SOLUTION["stress"], axis=0) / dt[:, None]
    solution_diff = np.log10(np.abs(solution_diff))

    plt.figure(figsize=(12, 9))
    plt.subplot(3, 2, 1)
    plt.plot(t_diff / SPY, solution_diff[:, 0::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}\dot{\sigma}_{xx}$ (Pa)")

    plt.subplot(3, 2, 2)
    plt.plot(solution_diff[:, 0::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}\dot{\sigma}_{xx}$ (Pa)")

    plt.subplot(3, 2, 3)
    plt.plot(t_diff / SPY, solution_diff[:, 1::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}\dot{\sigma}_{yy}$ (Pa)")

    plt.subplot(3, 2, 4)
    plt.plot(solution_diff[:, 1::3], linewidth=0.5)
    plt.ylabel(r"$\log_{10}\dot{\sigma}_{yy}$ (Pa)")

    plt.subplot(3, 2, 5)
    plt.plot(t_diff / SPY, solution_diff[:, 2::3], linewidth=0.5)
    plt.xlabel("time (years)")
    plt.ylabel(r"$\log_{10}\dot{\sigma}_{xy}$ (Pa)")

    plt.subplot(3, 2, 6)
    plt.plot(solution_diff[:, 2::3], linewidth=0.5)
    plt.xlabel("steps")
    plt.ylabel(r"$\log_{10}\dot{\sigma}_{xy}$ (Pa)")
    plt.suptitle("Stress rate changes")
    plt.show(block=False)


plot_stress_time_series_velocity()


def plot_invariants_time_series_velocity():
    """ Plot time integrated time series for each node """
    t_diff = SOLUTION["t"][0:-1]  # WHAAA???
    dt = np.diff(SOLUTION["t"])
    solution_diff = np.diff(SOLUTION["stress"], axis=0) / dt[:, None]

    # Calculate a few stress invariants
    I1 = solution_diff[:, 0::3] + solution_diff[:, 1::3]
    I2 = (
        solution_diff[:, 0::3] * solution_diff[:, 1::3] - solution_diff[:, 2::3] ** 2
    )  # 2nd invariant
    J2 = (I1 ** 2) / 3.0 - I2  # 2nd invariant (deviatoric)
    I1 = np.log(np.abs(I1))
    I2 = np.log(np.abs(I2))
    J2 = np.log(np.abs(J2))

    plt.figure(figsize=(12, 9))
    plt.subplot(3, 2, 1)
    plt.plot(t_diff / SPY, I1, linewidth=0.5)
    plt.ylabel(r"$\log_{10} \, |\dot{\mathrm{I}}_1|$ (Pa)")

    plt.subplot(3, 2, 2)
    plt.plot(I1, linewidth=0.5)
    plt.ylabel(r"$\log_{10}|\dot{\mathrm{I}}_1|$ (Pa)")

    plt.subplot(3, 2, 3)
    plt.plot(t_diff / SPY, I2, linewidth=0.5)
    plt.ylabel(r"$\log_{10} \, |\dot{\mathrm{I}}_2|$ (Pa$^2$)")

    plt.subplot(3, 2, 4)
    plt.plot(I2, linewidth=0.5)
    plt.ylabel(r"$\log_{10} \, |\dot{\mathrm{I}}_2|$ (Pa$^2$)")

    plt.subplot(3, 2, 5)
    plt.plot(t_diff / SPY, J2, linewidth=0.5)
    plt.xlabel("time (years)")
    plt.ylabel(r"$\log_{10} \, |\dot{\mathrm{J}}_2|$ (Pa$^2$)")

    plt.subplot(3, 2, 6)
    plt.plot(J2, linewidth=0.5)
    plt.xlabel("steps")
    plt.ylabel("$\log_{10} \, |\mathrm{\dot{J}}_2|$ (Pa$^2$)")
    plt.suptitle(r"Some stress invariant rate changes")
    plt.show(block=False)


plot_invariants_time_series_velocity()


def plot_slip_profile():
    """ Plot time integrated time series for each node """
    plot_times = np.floor(np.linspace(0, SOLUTION["y"].shape[0] - 1, 50)).astype(int)
    plot_x = np.array([_["x_integration_points"] for _ in ELEMENTS_FAULT]).flatten()
    plt.figure(figsize=(6, 9))
    y_labels = ["$u_x$ (m)", "$u_y$ (m)", "state"]
    for i, y_label in enumerate(y_labels):
        plt.subplot(3, 1, i + 1)
        for j in range(plot_times.size):
            plt.plot(plot_x, SOLUTION["y"][plot_times[j], i::3], linewidth=0.5)
        plt.ylabel(y_label)
    plt.xlabel("$x$ (m)")
    plt.show(block=False)


plot_slip_profile()


def plot_volume(time_idx, count):
    # Observation points for internal evaluation and visualization
    n_pts = 100
    x_plot = np.linspace(-25e3, 25e3, n_pts)
    y_plot = np.linspace(-25e3, 25e3, n_pts)
    x_plot, y_plot = np.meshgrid(x_plot, y_plot)
    x_plot = x_plot.flatten()
    y_plot = y_plot.flatten()
    obs_pts = np.array([x_plot, y_plot]).T.copy()

    def common_plot_elements():
        # Plot fault trace
        for element in ELEMENTS_FAULT:
            plt.plot(
                [element["x1"], element["x2"]],
                [element["y1"], element["y2"]],
                "-k",
                linewidth=1.0,
                zorder=50,
            )
        x_lim = np.array([x_plot.min(), x_plot.max()])
        y_lim = np.array([y_plot.min(), y_plot.max()])
        plt.xticks([x_lim[0], 0, x_lim[1]])
        plt.yticks([y_lim[0], 0, y_lim[1]])
        plt.gca().set_aspect("equal")
        plt.ylabel("$y$ (m)")

    # Internal evaluation for fault
    fault_slip = np.empty(2 * N_NODES)
    fault_slip[0::2] = SOLUTION["y"][time_idx, 0::3]
    fault_slip[1::2] = SOLUTION["y"][time_idx, 1::3]
    displacement_from_fault, stress_from_fault = bem2d.integrate(
        obs_pts, ELEMENTS_FAULT, PARAMETERS["mu"], PARAMETERS["nu"], "slip", fault_slip
    )

    fault_slip_prev = np.empty(2 * N_NODES)
    fault_slip_prev[0::2] = SOLUTION["y"][time_idx - 1, 0::3]
    fault_slip_prev[1::2] = SOLUTION["y"][time_idx - 1, 1::3]
    displacement_from_fault_prev, stress_from_fault_prev = bem2d.integrate(
        obs_pts,
        ELEMENTS_FAULT,
        PARAMETERS["mu"],
        PARAMETERS["nu"],
        "slip",
        fault_slip_prev,
    )

    # Backward difference to get velocities
    dt = SOLUTION["t"][time_idx] - SOLUTION["t"][time_idx - 1]
    displacement_from_fault = (
        displacement_from_fault_prev - displacement_from_fault
    ) / dt
    stress_from_fault_ = (stress_from_fault_prev - stress_from_fault) / dt

    ux_plot = displacement_from_fault[0, :]
    uy_plot = displacement_from_fault[1, :]
    u_plot_field = np.sqrt(ux_plot ** 2 + uy_plot ** 2)  # displacement magnitude
    u_plot_field = np.log10(u_plot_field)

    sxx_plot = stress_from_fault_[0, :]
    syy_plot = stress_from_fault_[1, :]
    sxy_plot = stress_from_fault_[2, :]
    I1 = sxx_plot + syy_plot  # 1st invariant
    I2 = sxx_plot * syy_plot - sxy_plot ** 2  # 2nd invariant
    J2 = (I1 ** 2) / 3.0 - I2  # 2nd invariant (deviatoric)
    s_plot_field = np.log10(np.abs(J2))

    contour_vec_1 = np.arange(-10, 2, 0.25)
    plt.figure(figsize=(8, 8))
    plt.subplot(2, 1, 1)
    plt.contourf(
        x_plot.reshape(n_pts, n_pts),
        y_plot.reshape(n_pts, n_pts),
        u_plot_field.reshape(n_pts, n_pts),
        # n_contours,
        contour_vec_1,
        extend="both",
        cmap=plt.get_cmap("plasma"),
    )
    plt.colorbar(fraction=0.046, pad=0.04, extend="both", label="$||v_i||$ (m)")
    plt.contour(
        x_plot.reshape(n_pts, n_pts),
        y_plot.reshape(n_pts, n_pts),
        u_plot_field.reshape(n_pts, n_pts),
        # n_contours,
        contour_vec_1,
        extend="both",
        linewidths=0.25,
        colors="k",
    )
    common_plot_elements()
    plt.gca().set_aspect("equal")

    n_contours = 20
    # contour_vec_2 = np.arange(-10, 11, 1)
    contour_vec_2 = np.arange(-20, 41, 1)

    plt.subplot(2, 1, 2)
    plt.contourf(
        x_plot.reshape(n_pts, n_pts),
        y_plot.reshape(n_pts, n_pts),
        s_plot_field.reshape(n_pts, n_pts),
        contour_vec_2,
        # n_contours,
        extend="both",
        cmap=plt.get_cmap("hot_r"),
    )
    plt.colorbar(
        fraction=0.046,
        pad=0.04,
        extend="both",
        label="$\log_{10} \, |\dot{\mathrm{J}}_2|$ (Pa$^2$)",
    )
    plt.contour(
        x_plot.reshape(n_pts, n_pts),
        y_plot.reshape(n_pts, n_pts),
        s_plot_field.reshape(n_pts, n_pts),
        contour_vec_2,
        # n_contours,
        extend="both",
        linewidths=0.25,
        colors="k",
    )
    common_plot_elements()
    plt.gca().set_aspect("equal")
    plt.xlabel("$x$ (m)")
    current_time = SOLUTION["t"][time_idx] / SPY
    plt.suptitle("t = " + f"{current_time:.4f}" + " (years)")
    # f"{SOLVER_TIME:.2f}" + " (seconds)")

    plt.show(block=False)
    print("Writing: " + OUTDIR + str(count).zfill(7) + ".png")
    plt.savefig(OUTDIR + str(count).zfill(7) + ".png")
    plt.close("all")


# n_frames = 50
# plot_idx = np.floor(np.linspace(1, SOLUTION["y"].shape[0] - 1, n_frames)).astype(int)
# plot_idx = np.floor(np.linspace(1, SOLUTION["y"].shape[0] - 1, n_frames)).astype(int)
# plot_idx = np.arange(1, SOLUTION["t"].size, 2)
# plot_idx = np.arange(10000, 29000, 100)
# for i in range(0, plot_idx.size):
#     plot_volume(plot_idx[i], i)

# # Convert .pngs to .mp4
# # TODO: make Quicktime compatible
# mp4_name = OUTDIR + str(time.time()) + ".mp4"
# png_name = OUTDIR + "%07d.png"
# convert_command = (
#     "ffmpeg -framerate 10 -i "
#     + png_name
#     + " -c:v libx264 -r 30 -y -v 32 -loglevel panic "
#     + mp4_name
# )
# os.system(convert_command)
