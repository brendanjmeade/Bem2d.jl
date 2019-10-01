using PyPlot

function fig_elementgeometry()
    fontsize = 10

    srcx = [-2, 2]
    srcy = [0, 0]
    obsx = -1.5
    obsy = -1.5

    angle = -35
    dx = 3.0
    dy = 2.3
    rotmat = [cosd(angle) -sind(angle); sind(angle) cosd(angle)]
    obsxglobal, obsyglobal = rotmat * [obsx ; obsy] + [dx ; dy]
    tx1, ty1 = rotmat * [srcx[1] ; srcy[1]] + [dx ; dy]
    tx2, ty2 = rotmat * [srcx[2] ; srcy[2]] + [dx ; dy]
    srcxglobal = [tx1, tx2]
    srcyglobal = [ty1, ty2]

    close("all")
    fig = figure(figsize=(5, 5))
    text(5.00, -0.30, L"$x$", fontsize=fontsize, ha="center", va="center")
    text(-0.30, 4.00, L"$y$", fontsize=fontsize, ha="center", va="center")

    plot(srcx, srcy, "-r", linewidth = 7, zorder=30, solid_capstyle="round")
    plot(srcx, srcy, ".w", linewidth = 0, zorder=30, solid_capstyle="round")
    text(-2.00, -0.5, L"$(x^{src}_1=-a, y^{src}_1=0)$", color="red", fontsize=fontsize, ha="center", va="center")
    text(2.00, -0.5, L"$(x^{src}_2=a, y^{src}_2=0)$", color="red", fontsize=fontsize, ha="center", va="center")
    plot(obsx, obsy, "ro", linewidth = 4, zorder=30)
    text(obsx, obsy-0.4, L"$(x^{obs}, y^{obs})$", color="red", fontsize=fontsize, ha="center", va="center")

    plot(srcxglobal, srcyglobal, "-b", linewidth = 7, zorder=30, solid_capstyle="round")
    plot(srcxglobal, srcyglobal, ".w", linewidth = 0, zorder=30, solid_capstyle="round")
    text(srcxglobal[1], srcyglobal[1] + 0.50, L"$(x^{src}_1, y^{src}_1)$", color="blue", fontsize=fontsize, ha="center", va="center")
    text(srcxglobal[2], srcyglobal[2] - 0.50, L"$(x^{src}_2, y^{src}_2)$", color="blue", fontsize=fontsize, ha="center", va="center")
    plot(obsxglobal, obsyglobal, "bo", linewidth = 4, zorder=30)
    text(obsxglobal, obsyglobal-0.5, L"$(\hat{x}^{obs}, \hat{y}^{obs})$", color="blue", fontsize=fontsize, ha="center", va="center")

    ax = gca()
    ax.set_xlim([-5, 5])
    ax.set_ylim([-4, 4])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(false) # Hide the top edge of the axis
    ax.spines["right"].set_visible(false) # Hide the right edge of the axis
    ax.spines["left"].set_position("center") # Move the right axis to the center
    ax.spines["bottom"].set_position("center") # Most the bottom axis to the center
    gca().set_aspect("equal")
    show()
end
fig_elementgeometry()
