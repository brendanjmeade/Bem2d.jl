using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using TriangleMesh
using Bem2d


"""
   kelvinUD(x, y, xoffset, yoffset, fx, fy, mu, nu)

   Analytic point source Kelvin solution from Crouch and Starfield
   section 4.2.
"""
function kelvinUD(x, y, xoffset, yoffset, fx, fy, mu, nu)
    x = x .- xoffset
    y = y .- yoffset
    Ukelvin = zeros(length(x), 2)
    Skelvin = zeros(length(x), 3)
    C = 1/(4*pi*(1-nu))
    r = sqrt.(x.^2 .+ y.^2)
    g = -C .* log.(r)
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
    # x, y = obsgrid(1e-4, 1e-4, 1-1e-4, 1-1e-4, npts)
    x, y = obsgrid(-1, -1, 1, 1, npts)

    fx = 0.0
    fy = 9.81 * 2700.0

    #! BEM geometry
    els = Elements(Int(1e5))
    # x1, y1, x2, y2 = discretizedline(0, 0, 1, 0, nels)
    x1, y1, x2, y2 = discretizedline(-1, -1, 1, -1, nels)
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

    # x1, y1, x2, y2 = discretizedline(1, 0, 1, 1, nels)
    x1, y1, x2, y2 = discretizedline(1, -1, 1, 1, nels)
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

    # x1, y1, x2, y2 = discretizedline(1, 1, 0, 1, nels)
    x1, y1, x2, y2 = discretizedline(1, 1, -1, 1, nels)
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

    # x1, y1, x2, y2 = discretizedline(0, 1, 0, 0, nels)
    x1, y1, x2, y2 = discretizedline(-1, 1, -1, -1, nels)
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
    mesh = trimesh(nodes, connectivity, 5)

    #! Loop over each triangle calculate area and contribution to u_eff
    #! For each element centroid calculate the effect of a body force
    ntri = size(mesh.cell)[2]
    Uboundarygravity = zeros(els.endidx, 2)
    for i in 1:ntri
        xtri = mesh.point[1, mesh.cell[:, i]]
        ytri = mesh.point[2, mesh.cell[:, i]]
        triarea = trianglearea(xtri[1], ytri[1], xtri[2], ytri[3], xtri[3], ytri[3])
        tricentroid = [mean(xtri) mean(ytri)]
        Ugi, _ = kelvinUD(els.xcenter[1:1:els.endidx],
                          els.ycenter[1:1:els.endidx],
                          tricentroid[1],
                          tricentroid[2],
                          fx,
                          fy,
                          mu,
                          nu)
        Uboundarygravity += triarea .* Ugi
    end
    
    #! What I need is a full forward evaluation at the bottom coordinates only
    Uboundarybottom, _ = constdispstress(slip2dispstress,
                                       els.xcenter[idx["bottom"]],
                                       els.ycenter[idx["bottom"]],
                                       els,
                                       collect(1:1:4*nels),
                                       Uboundarygravity[:, 1],
                                       Uboundarygravity[:, 2],
                                       mu,
                                       nu)

    #! Forward evaluation
    Uvolumegravity, Svolumegravity = constdispstress(slip2dispstress,
                                x,
                                y,
                                els,
                                collect(1:1:4*nels),
                                Uboundarygravity[:, 1],
                                Uboundarygravity[:, 2],
                                mu,
                                nu)
    Uvolumebottom, _ = constdispstress(slip2dispstress,
                                      x,
                                      y,
                                      els,
                                      idx["bottom"],
                                      Uboundarybottom[:, 1],
                                      Uboundarybottom[:, 2],
                                      mu,
                                      nu)
    # @show Uboundarybottom[:, 1]
    # @show asdf = reshape(Uvolumebottom[:, 1], npts, npts)

    # @infiltrate
    # return

    #! Plot meshes and fields
    figure(figsize=(20,20))
    nrows = 3
    ncols = 3
    fontsize = 20
    contour_levels = collect(LinRange(-2e-8, 10e-8, 30))
    contour_color = "white"
    contour_linewidth = 0.5
    markersize = 12

    subplot(nrows, ncols, 1)
    plotmesh(mesh)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("edge and area meshes", fontsize=fontsize)
    plot([x1, x2], [y1, y2], "-k", linewidth=2)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    subplot(nrows, ncols, 2)
    plot(Uboundarygravity[:, 1], "-r", label="x (Uboundarygravity)")
    plot(Uboundarygravity[:, 2], "-b", label="y (Uboundarygravity)")
    plot(Uboundarybottom[:, 1], "xr", label="x (Uboundarybottom)")
    plot(Uboundarybottom[:, 2], ".b", label="y (Uboundarybottom)")
    legend(fontsize=fontsize)
    title("boundary displacements", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)

    #! y-components only
    subplot(nrows, ncols, 3) # Plot slices
    xmat = reshape(x, npts, npts)
    ymat = reshape(y, npts, npts)
    uymat = reshape(Uvolumegravity[:, 2], npts, npts)
    uymatdispgbottom = reshape(Uvolumebottom[:, 2], npts, npts)
    plot(xmat[:, 1], uymat[:, 1], "-r", label="Uvolumegravity", linewidth=2.0)
    plot(xmat[:, 1], 2 .* uymatdispgbottom[:, 1], "r+", markersize=markersize, label="Uvolumebottom")
    plot(els.xcenter[idx["bottom"]], Uboundarybottom[:, 2], "bx", markersize=markersize, label="Uboundarybottom")
    gca().tick_params(labelsize=fontsize)
    legend(fontsize=fontsize)
    title("values at bottom of obs mat", fontsize=fontsize)

    subplot(nrows, ncols, 7)
    mat = Uvolumegravity[:, 2] #.- minimum(Uvolumegravity[:, 2]) 
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("Uvolumegravity", fontsize=fontsize)

    subplot(nrows, ncols, 8)
    mat = Uvolumebottom[:, 2] #.- minimum(Udispgbottom[:, 2])
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("Uvolumebottom", fontsize=fontsize)

    subplot(nrows, ncols, 9)
    mat1 = Uvolumegravity[:, 2] #.- minimum(Uvolumegravity[:, 2])
    mat2 = 2 .* Uvolumebottom[:, 2] #.- minimum(Uvolumebottom[:, 2])
    mat = mat1 .- mat2
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("uy", fontsize=fontsize)
    tight_layout()
    show()

    #! Classic 6-panel components only
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), Uvolumegravity, Svolumegravity, "BEM gravity only")
    subplot(2, 3, 1) # This overwrites a previous plot...ugly!
    plotmesh(mesh)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("edge and area meshes", fontsize=fontsize)
    plot([x1, x2], [y1, y2], "-k", linewidth=2)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    show()

end
gravitydisc()
