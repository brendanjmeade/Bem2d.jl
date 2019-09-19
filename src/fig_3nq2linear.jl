using Revise
using PyCall
using PyPlot

function ϕcoef(x, y, a)
    ∂ = zeros(length(x), 3)
    ∂[:, 1] = @. (x / a) * (9 * (x / a) / 8 - 3 / 4)
    ∂[:, 2] = @. (1 - 3 * (x / a) / 2) * (1 + 3 * (x / a) / 2)
    ∂[:, 3] = @. (x / a) * (9 * (x / a) / 8 + 3 / 4)
    return inv(∂) * y
end

function fig_3nq2linear()
    linewidth = 0.5
    fontsize = 6
    n = 1000
    a = 1
    x = LinRange(-a, a, n)

    # ϕ shape functions and coefficients
    ϕ1 = @. (x / a) * (9 * (x / a) / 8 - 3 / 4)
    ϕ2 = @. (1 - 3 * (x / a) / 2) * (1 + 3 * (x / a) / 2)
    ϕ3 = @. (x / a) * (9 * (x / a) / 8 + 3 / 4)
    xvals = [-1.0 ; 0.0 ; 1.0]
    yvals = @. [-1.0 ; 0.5 ; 2.0] + 1.25 # linear slip
    coef = ϕcoef(xvals, yvals, a)

    close("all")
    figure(figsize=(3, 3))
    plot(x, coef[1] * ϕ1, "--k", color="blue", linewidth=linewidth)
    plot(x, coef[2] * ϕ2, "--k", color="blue", linewidth=linewidth)
    plot(x, coef[3] * ϕ3, "--k", color="blue", linewidth=linewidth)
    plot(x, coef[1] * ϕ1 + coef[2] * ϕ2 + coef[3] * ϕ3, "-r", linewidth=linewidth)
    xticks([-1, 0, 1])
    yticks([-4, 0, 4])
    xlim([-1, 1])
    ylim([-4, 4])
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$u$ (m)", fontsize=fontsize)
    gca().tick_params(labelsize = fontsize)
    show()
end
fig_3nq2linear()
