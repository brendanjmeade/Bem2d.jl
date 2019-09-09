using Revise
using PyCall
using PyPlot
using Bem2d

function ex_constquadpartials()
    # Material and geometric constants
    μ = 3e10
    ν = 0.25
    nels = 10
    els = Elements(Int(1e5))
    L = 10000
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nels)
    els.x1[els.endidx + 1:els.endidx + nels] = x1
    els.y1[els.endidx + 1:els.endidx + nels] = y1
    els.x2[els.endidx + 1:els.endidx + nels] = x2
    els.y2[els.endidx + 1:els.endidx + nels] = y2
    els.name[els.endidx + 1:els.endidx + nels] .= "fault"
    standardize_elements!(els)

    srcidx = findall(x->x == "fault", els.name)
    obsidx = findall(x->x == "fault", els.name)
    ∂uconst, ∂σconst, ∂tconst = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)
    # ∂uquad, ∂σquad, ∂tquad = ∂quaduσ(slip2uσ, els, srcidx, obsidx, μ, ν)

    # Evaluation points and slip
    # xeval = [x["xnodes"] for x in elements][:]
    # yeval = [x["ynodes"] for x in elements][:]
    slipquad = zeros(6 * nels)
    slipconst = zeros(2 * nels)

    # Strike-slip
    slipquad[1:2:end] .= 1  # constant x-slip global
    slipconst[1:2:end] .= 1  # constant x-slip global
    # Tensile-slip
    # slipquad[2:2:end] .= 1  # constant x-slip global
    # slipconst[2:2:end] .= 1  # constant x-slip global

    # Predict on fault displacements, stresses, and tractions
    uconst = ∂uconst * slipconst
    σconst = ∂σconst * slipconst
    tconst = ∂tconst * slipconst
    # uquad = ∂uquad * slipquad
    # σquad = ∂σquad * slipquad
    # tquad = ∂tquad * slipquad

    # # Plot geometry of elements
    # plt.figure(figsize=(12, 8))
    # plt.subplot(2, 2, 1)
    # for element in elements:
    #     plt.plot(
    #         [element["x1"], element["x2"]],
    #         [element["y1"], element["y2"]],
    #         "-k",
    #         color="r",
    #         linewidth=0.5,
    #     )
    #     plt.plot(
    #         [element["x1"], element["x2"]],
    #         [element["y1"], element["y2"]],
    #         "r.",
    #         markersize=1,
    #         linewidth=0.5,
    #     )
    # for i, element in enumerate(elements):
    #     plt.text(
    #         element["x_center"],
    #         element["y_center"],
    #         str(i),
    #         horizontalalignment="center",
    #         verticalalignment="center",
    #         fontsize=8,
    #     )
    # plt.xlabel("x (m)")
    # plt.ylabel("y (m)")
    # plt.title("element geometry")

    # # Plot predicted displacements
    # plt.subplot(2, 2, 2)
    # plt.plot(
    #     x_eval[1::3],
    #     displacement_quadratic[2::6],
    #     "+r",
    #     markeredgewidth=0.5,
    #     label="u_x quadratic",
    # )
    # plt.plot(
    #     x_eval[1::3],
    #     displacement_quadratic[3::6],
    #     "+b",
    #     markeredgewidth=0.5,
    #     label="u_y quadratic",
    # )
    # plt.plot(
    #     x_eval[1::3],
    #     displacement_constant[0::2],
    #     "or",
    #     markerfacecolor="none",
    #     markeredgewidth=0.5,
    #     label="u_x constant",
    # )
    # plt.plot(
    #     x_eval[1::3],
    #     displacement_constant[1::2],
    #     "ob",
    #     markerfacecolor="none",
    #     markeredgewidth=0.5,
    #     label="u_y constant",
    # )
    # plt.legend(loc="upper right")
    # plt.xlabel("x (m)")
    # plt.ylabel("displacements (m)")
    # plt.title("displacements")

    # plt.subplot(2, 2, 3)
    # plt.plot(x_eval, slip_quadratic[0::2], "+k", label="quadratic", markeredgewidth=0.5)
    # plt.plot(
    #     x_eval[1::3],
    #     slip_constant[0::2],
    #     "ok",
    #     markerfacecolor="none",
    #     label="constant",
    #     markeredgewidth=0.5,
    # )
    # plt.legend()
    # plt.xlabel("x (m)")
    # plt.ylabel("input slip")
    # plt.title("element slip")

    # plt.subplot(2, 2, 4)
    # plt.plot(
    #     x_eval, stress_quadratic[0::3], "+r", label="s_xx quadratic", markeredgewidth=0.5
    # )
    # plt.plot(
    #     x_eval, stress_quadratic[1::3], "+b", label="s_yy quadratic", markeredgewidth=0.5
    # )
    # plt.plot(
    #     x_eval, stress_quadratic[2::3], "+k", label="s_xy quadratic", markeredgewidth=0.5
    # )

    # plt.plot(
    #     x_eval[1::3],
    #     stress_constant[0::3],
    #     "or",
    #     markerfacecolor="none",
    #     markeredgewidth=0.5,
    #     label="s_xx constant",
    # )
    # plt.plot(
    #     x_eval[1::3],
    #     stress_constant[1::3],
    #     "ob",
    #     markerfacecolor="none",
    #     markeredgewidth=0.5,
    #     label="s_yy constant",
    # )
    # plt.plot(
    #     x_eval[1::3],
    #     stress_constant[2::3],
    #     "ok",
    #     markerfacecolor="none",
    #     markeredgewidth=0.5,
    #     label="s_xy constant",
    # )
    # plt.legend()
    # plt.xlabel("x (m)")
    # plt.ylabel("stresses (Pa)")
    # plt.title("stresses")
    # plt.show(block=False)
end
ex_constquadpartials()
