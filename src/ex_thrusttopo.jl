using Revise
using Plots
using Bem2d
pyplot()


function ex_thrusttopo()
    # TOPO_SIGN_FLIP = True
    mu = 30e9
    nu = 0.25
    elements = Elements()

    # Observation points for internal evaluation and visualization
    npts = 50
    xobs = LinRange(-10e3, 10e3, npts)
    yobs = LinRange(-5e3, 5e3, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # Topographic free surface
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, 20)
    y1 = -1e3 * @.atan(x1 / 1e3)
    y2 = -1e3 * @.atan(x2 / 1e3)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "freesurface"
    end
    standardize_elements!(elements)

    # Curved fault
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, 10)
    y1 = 3e3 * @.atan(x1 / 1e3)
    y2 = 3e3 * @.atan(x2 / 1e3)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "fault"
    end
    standardize_elements!(elements)

    # Partial derivatves
    srcidx = findall(x -> x == "fault", elements.name)
    obsidx = findall(x -> x == "freesurface", elements.name)
    d1, s1, t1 = partials_constslip(elements, srcidx, obsidx, mu, nu)
    srcidx = findall(x -> x == "freesurface", elements.name)
    obsidx = findall(x -> x == "freesurface", elements.name)
    d2, s2, t2 = partials_constslip(elements, srcidx, obsidx, mu, nu)

    # Remove and separate BCs, local to global transform here?
    nfaultelements = length(findall(x -> x == "fault", elements.name))
    faultslip = zeros(2 * nfaultelements)
    for i in 1:nfaultelements # There has to be a matlab style way to do this without a loop
       faultslip[2 * (i - 1) + 1] = 1.0  # TODO: This is in global not local!!!
    end

    # Solve the BEM problem
    disp_freesurface = inv(t2) * t1 * faultslip

    # Fault in full space
    dispfault = zeros(length(xobs), 2)
    stressfault = zeros(length(xobs), 3)
    faultidx = findall(x -> x == "fault", elements.name)
    for i in faultidx[1]:faultidx[end]
        disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[i],
            mu, nu, 1, 0, elements.xcenter[i], elements.ycenter[i],
            elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
        dispfault += disp
        stressfault += stress
    end
    plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        dispfault, stressfault, "fault")

    # Free surface in full space
    # TODO: Need to adopt partials method of indexing into indices
    # dispfault = zeros(length(xobs), 2)
    # stressfault = zeros(length(xobs), 3)
    # faultidx = findall(x -> x == "fault", elements.name)
    # for i in faultidx[1]:faultidx[end]
    #     disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[i],
    #         mu, nu, 1, 0, elements.xcenter[i], elements.ycenter[i],
    #         elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
    #     dispfault += disp
    #     stressfault += stress
    # end
    # plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
    #     dispfault, stressfault, "fault")


    # displacement_from_topo, stress_from_topo = bem2d.integrate(
    #     obs_pts, elements_surface, mu, nu, "slip", displacement_free_surface
    # )
    #
    # # Pretty of displacements and stresses
    # def common_plot_elements():
    #     # Create a white fill over portion of the figure above the free surface
    #     x_surface = np.unique(
    #         [[_["x1"] for _ in elements_surface], [_["x2"] for _ in elements_surface]]
    #     )
    #     x_fill = np.append(x_surface, [10e3, -10e3, -10e3])
    #     y_surface = np.unique(
    #         [[_["y1"] for _ in elements_surface], [_["y2"] for _ in elements_surface]]
    #     )
    #     y_surface = np.flip(y_surface, 0)
    #     y_fill = np.append(y_surface, [5e3, 5e3, np.min(y_surface)])
    #     plt.fill(x_fill, y_fill, "w", zorder=30)
    #
    #     for element in elements_fault + elements_surface:
    #         plt.plot(
    #             [element["x1"], element["x2"]],
    #             [element["y1"], element["y2"]],
    #             "-k",
    #             linewidth=1.0,
    #             zorder=50,
    #         )
    #
    #     x_lim = np.array([x_plot.min(), x_plot.max()])
    #     y_lim = np.array([y_plot.min(), y_plot.max()])
    #     plt.xticks([x_lim[0], 0, x_lim[1]])
    #     plt.yticks([y_lim[0], 0, y_lim[1]])
    #     plt.gca().set_aspect("equal")
    #     plt.xlabel("$x$ (m)")
    #     plt.ylabel("$y$ (m)")
    #
    #
    # if TOPO_SIGN_FLIP:
    #     displacement_from_topo *= -1
    #     stress_from_topo *= -1
    #
    #
    # ux_plot = (displacement_from_topo + displacement_from_fault)[0, :]
    # uy_plot = (displacement_from_topo + displacement_from_fault)[1, :]
    # u_plot_field = np.sqrt(ux_plot ** 2 + uy_plot ** 2)  # displacement magnitude
    #
    # sxx_plot = (stress_from_topo + stress_from_fault)[0, :]
    # syy_plot = (stress_from_topo + stress_from_fault)[1, :]
    # sxy_plot = (stress_from_topo + stress_from_fault)[2, :]
    # I1 = sxx_plot + syy_plot  # 1st invariant
    # I2 = sxx_plot * syy_plot - sxy_plot ** 2  # 2nd invariant
    # J2 = (I1 ** 2) / 3.0 - I2  # 2nd invariant (deviatoric)
    # s_plot_field = np.log10(np.abs(J2))
    #
    #
    # n_contours = 5
    # plt.figure(figsize=(6, 8))
    # plt.subplot(2, 1, 1)
    # plt.contourf(
    #     x_plot.reshape(n_pts, n_pts),
    #     y_plot.reshape(n_pts, n_pts),
    #     u_plot_field.reshape(n_pts, n_pts),
    #     n_contours,
    #     cmap=plt.get_cmap("plasma"),
    # )
    # plt.colorbar(fraction=0.046, pad=0.04, extend="both", label="$||u_i||$ (m)")
    # plt.contour(
    #     x_plot.reshape(n_pts, n_pts),
    #     y_plot.reshape(n_pts, n_pts),
    #     u_plot_field.reshape(n_pts, n_pts),
    #     n_contours,
    #     linewidths=0.25,
    #     colors="k",
    # )
    # common_plot_elements()
    # plt.title("displacement magnitude")
    #
    # plt.subplot(2, 1, 2)
    # plt.contourf(
    #     x_plot.reshape(n_pts, n_pts),
    #     y_plot.reshape(n_pts, n_pts),
    #     s_plot_field.reshape(n_pts, n_pts),
    #     n_contours,
    #     cmap=plt.get_cmap("hot_r"),
    # )
    # plt.colorbar(
    #     fraction=0.046, pad=0.04, extend="both", label="$log_{10}|\mathrm{J}_2|$ (Pa$^2$)"
    # )
    # plt.contour(
    #     x_plot.reshape(n_pts, n_pts),
    #     y_plot.reshape(n_pts, n_pts),
    #     s_plot_field.reshape(n_pts, n_pts),
    #     n_contours,
    #     linewidths=0.25,
    #     colors="k",
    # )
    # common_plot_elements()
    # plt.title("second stress invariant (deviatoric)")
    # plt.show(block=False)
    #
    #
    # bem2d.plot_fields(
    #     elements_surface + elements_fault,
    #     x_plot.reshape(n_pts, n_pts),
    #     y_plot.reshape(n_pts, n_pts),
    #     displacement_from_fault,
    #     stress_from_fault,
    #     "Fault",
    # )
    #
    # bem2d.plot_fields(
    #     elements_surface + elements_fault,
    #     x_plot.reshape(n_pts, n_pts),
    #     y_plot.reshape(n_pts, n_pts),
    #     displacement_from_topo,
    #     stress_from_topo,
    #     "Topography",
    # )
    #
    # bem2d.plot_fields(
    #     elements_surface + elements_fault,
    #     x_plot.reshape(n_pts, n_pts),
    #     y_plot.reshape(n_pts, n_pts),
    #     displacement_from_topo + displacement_from_fault,
    #     stress_from_topo + stress_from_fault,
    #     "Topography + fault",
    # )
end
ex_thrusttopo()
