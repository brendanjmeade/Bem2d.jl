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
        # mesh = refine(mesh, divide_cell_into=refinefactor, keep_edges=true)
        mesh = refine(mesh, divide_cell_into=refinefactor, keep_edges=false)
    end
    return mesh
end


"""
    plotmesh(mesh)

    Plot mesh with shading, mesh edges, and centroids
"""
function plotmesh(mesh)
    for i in 1:size(mesh.cell)[2]
        x = mesh.point[1, mesh.cell[:, i]]
        y = mesh.point[2, mesh.cell[:, i]]
        fill(x, y, "lightgray", edgecolor="k", linewidth=0.5)
    end
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
    nels = 20
    npts = 50
    # x, y = obsgrid(0, 0, 1, 1, npts)
    x, y = obsgrid(1e-4, 1e-4, 1-1e-4, 1-1e-4, npts)
    fx = 0.0
    fy = 9.81 * 2700.0

    #! BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(0, 0, 1, 0, nels)
    xbottom = x1
    ybottom = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "bottom"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(1, 0, 1, 1, nels)
    xright = x1
    yright = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "right"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(1, 1, 0, 1, nels)
    xtop = x1
    ytop = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "top"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(0, 1, 0, 0, nels)
    xleft = x1
    yleft = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "left"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Discretize disc volume into cells
    xcoords = [xbottom;xright;xtop;xleft] 
    ycoords = [ybottom;yright;ytop;yleft]
    nodes = [xcoords ycoords]
    connectivity = Int.(zeros(size(nodes)))
    connectivity[:, 1] = collect(1:1:length(xcoords))
    connectivity[1:end-1, 2] = collect(2:1:length(xcoords))
    connectivity[end, 2] = 1
    mesh = trimesh(nodes, connectivity, 1)

    #! Loop over each triangle calculate area and contribution to u_eff
    #! For each element centroid calculate the effect of a body force
    ntri = size(mesh.cell)[2]
    Ugmesh = zeros(els.endidx, 2)
    for i in 1:ntri
        xtri = mesh.point[1, mesh.cell[:, i]]
        ytri = mesh.point[2, mesh.cell[:, i]]
        triarea = trianglearea(xtri[1], ytri[1], xtri[2], ytri[3], xtri[3], ytri[3])
        tricentroid = [mean(xtri) mean(ytri)]
        Ugi, _ = kelvinUD(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], tricentroid[1], tricentroid[2], fx, fy, mu, nu)
        Ugmesh += @. triarea * Ugi
    end
    
    #! What I need is a full forward evaluation at the bottom coordinates only
    Udispgbottom2, _ = constdispstress(slip2dispstress,
                                       els.xcenter[idx["bottom"]],
                                       els.ycenter[idx["bottom"]],
                                       els,
                                       collect(1:1:4*nels),
                                       Ugmesh[:, 1],
                                       Ugmesh[:, 2],
                                       mu, nu)


    #! Forward evaluation
    Udispg, _ = constdispstress(slip2dispstress, x, y, els,
                                collect(1:1:4*nels),
                                Ugmesh[:, 1],
                                Ugmesh[:, 2],
                                mu, nu)
    Udispgbottom, _ = constdispstress(slip2dispstress, x, y, els,
                                      idx["bottom"],
                                      Udispgbottom2[:, 1],
                                      Udispgbottom2[:, 2],
                                      mu, nu)
    
    #! Plot meshes and fields
    figure(figsize=(20,20))
    nrows = 3
    ncols = 3
    fontsize = 20
    contour_levels = 20
    contour_color = "white"
    contour_linewidth = 1.0

    subplot(nrows, ncols, 1)
    plotmesh(mesh)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("edge and area meshes", fontsize=fontsize)
    plot([x1, x2], [y1, y2], "-k", linewidth=2)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 2)
    plot(Ugmesh[:, 1], "-r", label="ux gravity")
    plot(Ugmesh[:, 2], "-b", label="uy gravity")
    plot(Udispgbottom2[:, 1], "xr", label="ux correction")
    plot(Udispgbottom2[:, 2], ".b", label="uy correction")
    legend()
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 4)
    mat = Udispg[:, 1]
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("ux", fontsize=fontsize)

    subplot(nrows, ncols, 7)
    mat = Udispg[:, 2]
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    clim([0, 2e-8])
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("uy", fontsize=fontsize)

    subplot(nrows, ncols, 8)
    mat = Udispgbottom[:, 2]
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    clim([0, 2e-8])
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("uy", fontsize=fontsize)

    tight_layout()
    show()
end
gravitydisc()
