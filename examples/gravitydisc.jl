using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using TriangleMesh
using Bem2d


"""
    discretized_arc(θstart, θend, radius, n_pts)

Generate regularly spaced eleemnts along an curved arc.
"""
function discretized_arc(θstart, θend, radius, n_pts)
    # Create geometry of discretized arc
    θrange = collect(LinRange(θstart, θend, n_pts + 1))
    x = @. radius * cos(θrange)
    y = @. radius * sin(θrange)
    x1 = x[1:1:end-1]
    x2 = x[2:1:end]
    y1 = y[1:1:end-1]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end


"""
    circle_subplot(nrows, ncols, plotidx, x, y, mat, npts, R, theta0, title_string)

Plot field (displacement, stress) within a circular disk and style
"""
function circle_subplot(nrows, ncols, plotidx, x, y, mat, npts, R, theta0, title_string)
    fontsize = 20
    contour_levels = 200
    contour_color = "white"
    contour_linewidth = 0.5
    color_scale = 1

    subplot(nrows, ncols, plotidx)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title(title_string, fontsize=fontsize)

    # Draw entire circle and region of applied tractions
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, 360)
    plot([x1, x2], [y1, y2], "-k", linewidth=2)
    # x1, y1, x2, y2 = discretized_arc(-theta0, theta0, R, 50)
    # plot([x1, x2], [y1, y2], "-r", linewidth=2)
    # x1, y1, x2, y2 = discretized_arc(-theta0+deg2rad(180), theta0+deg2rad(180), R, 50)
    # plot([x1, x2], [y1, y2], "-r", linewidth=2)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
end


"""
   kelvinUD(x, y, xoffset, yoffset, fx, fy, mu, nu)

   Analytic point source Kelvin solution from Crouch and Starfield
   section 4.2.
"""
function kelvinUD(x, y, xoffset, yoffset, fx, fy, mu, nu)
    x = @. x - xoffset
    y = @. y - yoffset
    Ukelvin = zeros(length(x), 2)
    Skelvin = zeros(length(x), 3)
    C = 1/(4*pi*(1-nu))
    r = @. sqrt(x^2+y^2)
    g = @. -C * log(r)
    gx = @. -C * x/(x^2+y^2)
    gy = @. -C * y/(x^2+y^2)
    gxy = @. C * 2*x*y/(x^2+y^2)^2
    gxx = @. C * (x^2-y^2)/(x^2+y^2)^2
    gyy = -gxx
    Ukelvin[:, 1] = @. fx/(2*mu)*((3-4*nu)*g-x*gx) + fy/(2*mu)*(-y*gx)
    Ukelvin[:, 2] = @. fx/(2*mu)*(-x*gy) + fy/(2*mu)*((3-4*nu)*g-y*gy)
    Skelvin[:, 1] = @. fx*(2*(1-nu)*gx-x*gxx) + fy*(2*nu*gy-y*gxx)
    Skelvin[:, 2] = @. fx*(2*nu*gx-x*gyy) + fy*(2*(1-nu)*gy-y*gyy)
    Skelvin[:, 3] = @. fx*((1-2*nu)*gy-x*gxy) + fy*((1-2*nu)*gx-y*gxy)
    return Ukelvin, Skelvin
end


"""
    trianglearea(x1, y1, x2, y2, x3, y3)

    From https://math.stackexchange.com/questions/516219/finding-out-the-area-of-a-triangle-if-the-coordinates-of-the-three-vertices-are
"""
function trianglearea(x1, y1, x2, y2, x3, y3)
    return abs(0.5*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)))
end


"""
    trimesh(nodes, connectivity, refinefactor)

Calls triangle binary (https://www.cs.cmu.edu/~quake/triangle.html)
Based on example at: https://github.com/konsim83/TriangleMesh.jl/blob/master/docs/src/man/examples.md
"""
function trimesh(nodes, connectivity, refinefactor)
    nnodes = size(nodes)[1]
    poly = Polygon_pslg(nnodes, 1, 2, nnodes, 0)
    set_polygon_point!(poly, nodes)
    set_polygon_segment!(poly, connectivity)
    mesh = create_mesh(poly)
    if refinefactor > 1
        mesh = refine(mesh, divide_cell_into=refinefactor, keep_edges=true)
    end
    return mesh
end


#! Plot mesh
function plotmesh(mesh)
    figure()
    for i in 1:size(mesh.cell)[2]
        x = mesh.point[1, mesh.cell[:, i]]
        y = mesh.point[2, mesh.cell[:, i]]
        fill(x, y, "lightgray", edgecolor="k", linewidth=0.5)
        plot(mean(x), mean(y), "r.")
    end
    show()
    return nothing
end


"""
    gravitydisc()

Experiments with gravity body force.
"""
function gravitydisc()
    close("all")
    mu = 3e10
    nu = 0.25
    p = 1.0e-5 # Applied radial pressure over arc
    theta0 = deg2rad(15.0) # Arc length over which pressure is applied
    nels = 36
    R = nels / (2 * pi)
    npts = 300
    x, y = obsgrid(-R, -R, R, R, npts)
    r = @. sqrt(x^2 + y^2)

    #! BEM solution
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "circle"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Discretize disc volume into cells
    nodes = [x1 y1]
    connectivity = Int.(zeros(size(nodes)))
    connectivity[:, 1] = collect(1:1:length(x1))
    connectivity[1:end-1, 2] = collect(2:1:length(x1))
    connectivity[end, 2] = 1
    mesh = trimesh(nodes, connectivity, 1)
    plotmesh(mesh)

    #! Loop over each triangle calculate area and contribution to u_eff
    ntri = size(mesh.cell)[2]
    for i in 1:ntri
        @show i
        xtri = mesh.point[1, mesh.cell[:, i]]
        ytri = mesh.point[2, mesh.cell[:, i]]
        @show triarea = trianglearea(xtri[1], ytri[1], xtri[2], ytri[3], xtri[3], ytri[3])
        tricentroid = [mean(xtri) mean(ytri)]
        # Ug, _ = kelvinUD(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], 0, 0, fx, fy, mu, nu)
    end
    @infiltrate
    return

    #! Apply normal tractions everywhere and convert from radial to Cartesian
    xtracC = zeros(els.endidx)
    ytracC = zeros(els.endidx)
    thetaels = @. atan(els.ycenter[1:1:els.endidx], els.xcenter[1:1:els.endidx])
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtracC[i], ytracC[i] = els.rotmat[i, :, :] * normalTractions
    end

    #! Zero out the tractions on the area without contact
    deleteidx = findall(x -> (x>theta0 && x<deg2rad(180)-theta0), thetaels)
    xtracC[deleteidx] .= 0
    ytracC[deleteidx] .= 0
    deleteidx = findall(x -> (x<-theta0 && x>-deg2rad(180)+theta0), thetaels)
    xtracC[deleteidx] .= 0
    ytracC[deleteidx] .= 0


    #! For each element centroid calculate the effect of a body force
    fx = 0.0
    fy = 9.81 * 2700.0
    Ug, Sg = kelvinUD(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], 0, 0, fx, fy, mu, nu)
    Ug[9, 1] = 0.0 # Set a fixed point at bottom
    Ug[10, 1] = 0.0 # Set a fixed point at bottom
    Ug[9, 2] = 0.0 # Set a fixed point at bottom
    Ug[10, 2] = 0.0 # Set a fixed point at bottom

    figure()
    for i in 1:els.endidx
        text(els.xcenter[i], els.ycenter[i], string(i))
    end
    plot([els.x1[:] els.x2[:]], [els.y1[:] els.y2[:]], "-k")
    plot([x1, x2], [y1, y2], "-r", linewidth=2)
    show()

    figure()
    plot(Ug[:, 1], "r-", label="ux")
    plot(Ug[:, 2], "b-", label="uy")
    legend()

    #! Forward evaluation.  No BEM solve required for u_eff
    Udisp, Sdisp = constdispstress(slip2dispstress, x, y, els, idx["circle"], xtracC, ytracC, mu, nu)
    Udispg, Sdispg = constdispstress(slip2dispstress, x, y, els, idx["circle"], Ug[:, 1], Ug[:, 2], mu, nu)

    #! Isolate the values inside the circle
    nanidx = findall(x -> x > R, r)
    Udisp[nanidx, :] .= NaN
    Sdisp[nanidx, :] .= NaN
    Udispg[nanidx, :] .= NaN
    Sdispg[nanidx, :] .= NaN

    # figure()
    # Udispg[:, 2] = @. Udispg[:, 2] - minimum(Udispg[:, 2])
    # quiver(x, y, Udispg[:, 1],  Udispg[:, 2])

    #! 9-panel plot for quadratic only
    figure(figsize=(30,10))
    nrows = 1
    ncols = 5
    # circle_subplot(nrows, ncols, 1, x, y, Udisp[:, 1], npts, R, theta0, L"u_{x} \; \mathrm{(BEM)}")
    # circle_subplot(nrows, ncols, 2, x, y, Udisp[:, 2], npts, R, theta0, L"u_{y} \; \mathrm{(BEM)}")
    # circle_subplot(nrows, ncols, 3, x, y, Sdisp[:, 1], npts, R, theta0, L"\sigma_{xx} \; \mathrm{(BEM)}")
    # circle_subplot(nrows, ncols, 4, x, y, Sdisp[:, 2], npts, R, theta0, L"\sigma_{yy} \; \mathrm{(BEM)}")
    # circle_subplot(nrows, ncols, 5, x, y, Sdisp[:, 3], npts, R, theta0, L"\sigma_{xy} \; \mathrm{(BEM)}")
    circle_subplot(nrows, ncols, 1, x, y, Udispg[:, 1], npts, R, theta0, L"u_{x} \; \mathrm{(gravity BEM)}")
    circle_subplot(nrows, ncols, 2, x, y, Udispg[:, 2], npts, R, theta0, L"u_{y} \; \mathrm{(gravity BEM)}")
    circle_subplot(nrows, ncols, 3, x, y, Sdispg[:, 1], npts, R, theta0, L"\sigma_{xx} \; \mathrm{(gravity BEM)}")
    circle_subplot(nrows, ncols, 4, x, y, Sdispg[:, 2], npts, R, theta0, L"\sigma_{yy} \; \mathrm{(gravity BEM)}")
    circle_subplot(nrows, ncols, 5, x, y, Sdispg[:, 3], npts, R, theta0, L"\sigma_{xy} \; \mathrm{(gravity BEM)}")
    tight_layout()
    show()
end
gravitydisc()
