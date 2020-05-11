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
function kelvinUS(x, y, xoffset, yoffset, fx, fy, mu, nu)
    x = @. x - xoffset
    y = @. y - yoffset
    Ukelvin = zeros(length(x), 2)
    Skelvin = zeros(length(x), 3)
    C = 1/(4*pi*(1-nu))
    r = @. sqrt(x^2 + y^2)
    g = @. -C * log(r)
    gx = @. -C * x/(x^2+y^2)
    gy = @. -C * y/(x^2+y^2)
    gxy = @. C * 2*x*y/(x^2+y^2)^2
    gxx = @. C * (x^2-y^2)/(x^2+y^2)^2
    gyy = -gxx

    #! THERE IS AN ERROR IN HOW THIS WORKS FOR PARTIALS....WHAA?!?!?!
    Ukelvin[:, 1] = @. fx/(2*mu)*((3-4*nu)*g-x*gx) + fy/(2*mu)*(-y*gx)
    Ukelvin[:, 2] = @. fx/(2*mu)*(-x*gy) + fy/(2*mu)*((3-4*nu)*g-y*gy)
    Skelvin[:, 1] = @. fx*(2*(1-nu)*gx-x*gxx) + fy*(2*nu*gy-y*gxx)
    Skelvin[:, 2] = @. fx*(2*nu*gx-x*gyy) + fy*(2*(1-nu)*gy-y*gyy)
    Skelvin[:, 3] = @. fx*((1-2*nu)*gy-x*gxy) + fy*((1-2*nu)*gx-y*gxy)
    return Ukelvin, Skelvin
end


"""
    trianglearea(x1, y1, x2, y2, x3, y3)

    It's just half the absolute value of the determinant: 
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
    localsubplot(???)

    Common plotting commands
"""
function localsubplot(x, y, mat, titlestring, contourlabels)
    return nothing
end


"""
    PUTK(els, obsidx, srcidx, mu, nu)

Kernel relating body forces to boundary displacements and tractions.
Uses Kelvin kernels.
"""
function PUTK(els, obsidx, mesh, mu, nu)
    ntri = size(mesh.cell)[2]
    nobs = length(obsidx)
    Ku = zeros(2*nobs, 2*ntri)
    Kt = zeros(2*nobs, 2*ntri)
    _Ku = zeros(2, 2)
    _Ks = zeros(3, 2)
    _Kt = zeros(2, 2)
    for iobs in 1:nobs # Should I loop over obs or src?  src is easier to write.
        for isrc in 1:ntri
            # Triangle areas and centroids
            xtri = mesh.point[1, mesh.cell[:, isrc]]
            ytri = mesh.point[2, mesh.cell[:, isrc]]
            triarea = trianglearea(xtri[1], ytri[1], xtri[2], ytri[3], xtri[3], ytri[3])
            tricentroid = [mean(xtri) mean(ytri)]

            # Calculate displacement and stress kernels
            _Ku[:, 1], _Ks[:, 1] = kelvinUS([els.xcenter[obsidx[iobs]]], [els.ycenter[obsidx[iobs]]], tricentroid[1], tricentroid[2], 1, 0, mu, nu)
            _Ku[:, 2], _Ks[:, 2] = kelvinUS([els.xcenter[obsidx[iobs]]], [els.ycenter[obsidx[iobs]]], tricentroid[1], tricentroid[2], 0, 1, mu, nu)

            # Calculate traction kernels from stress kernels
            _Kt[:, 1] = stress2trac(_Ks[:, 1], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])
            _Kt[:, 2] = stress2trac(_Ks[:, 2], [els.xnormal[obsidx[iobs]] ; els.ynormal[obsidx[iobs]]])

            # Insert 2x2 matrices into system matrices
            Ku[2 * (iobs - 1) + 1:2 * (iobs - 1) + 2, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _Ku
            Kt[2 * (iobs - 1) + 1:2 * (iobs - 1) + 2, 2 * (isrc - 1) + 1:2 * (isrc - 1) + 2] = _Kt
        end
    end
    return Ku, Kt
end


"""
    gravitysquare()

Experiments with gravity body force.
"""
function gravitysquare()
    close("all")
    fontsize = 20
    mu = 3e10
    nu = 0.25
    nels = 20
    npts = 25
    L = 1e4
    x, y = obsgrid(-L+1, -L+1, L-1, L-1, npts)
    fx = 0.0
    fy = 9.81 * 2700.0

    #! BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-L, -L, L, -L, nels)
    xbottom = x1
    ybottom = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "B"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(L, -L, L, L, nels)
    xright = x1
    yright = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "R"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(L, L, -L, L, nels)
    xtop = x1
    ytop = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "T"
    end
    standardize_elements!(els)

    x1, y1, x2, y2 = discretizedline(-L, L, -L, -L, nels)
    xleft = x1
    yleft = y1
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "L"
    end
    standardize_elements!(els)
    idx = getidxdict(els)

    #! Triangulate into cells
    xcoords = [xbottom;xright;xtop;xleft] 
    ycoords = [ybottom;yright;ytop;yleft]
    nodes = [xcoords ycoords]
    connectivity = Int.(zeros(size(nodes)))
    connectivity[:, 1] = collect(1:1:length(xcoords))
    connectivity[1:end-1, 2] = collect(2:1:length(xcoords))
    connectivity[end, 2] = 1
    mesh = trimesh(nodes, connectivity, 1)
    ntri = size(mesh.cell)[2]
    
    #! Kernels matrices (boundary -> boundary)
    Bidx = idx["B"]
    RTLidx = [idx["R"] ; idx["T"]; idx["L"]]
    BRTLidx = collect(1:1:els.endidx)
    T_pB_qBRTL, H_pB_qBRTL = PUTC(slip2dispstress, els, Bidx, BRTLidx, mu, nu)
    T_pRTL_qBRTL, H_pRTL_qBRTL = PUTC(slip2dispstress, els, RTLidx, BRTLidx, mu, nu)

    #! Kernel matrices (volume -> boundary)
    Ku_pB_qv, Kt_pB_qv = PUTK(els, Bidx, mesh, mu, nu)
    Ku_pRTL_qv, Kt_pRTL_qv = PUTK(els, RTLidx, mesh, mu, nu)

    #! Assemble and solve BEM problem
    alpha = 7e-8
    TH = [T_pB_qBRTL ; alpha .* H_pRTL_qBRTL]
    K = [Ku_pB_qv ; alpha .* Kt_pRTL_qv]

    Tall = [T_pB_qBRTL ; alpha .* T_pRTL_qBRTL]
    Kall = [Ku_pB_qv ; alpha .* Ku_pRTL_qv]

    f = zeros(2*ntri)
    f[2:2:end] .= -9.8 * 2700
    Ueff = inv(TH) * K * f

    #! Does this Ueff actually give the correct boundary conditions?
    Ugravity = -K * f
    recoveredbcs = (TH * Ueff) + (-K * f)
    Uall = (Tall * Ueff) + (-Kall * f)
    Uboundaryonly = Tall * Ueff
    Ugravity = -Kall * f
    
    figure()
    plot(recoveredbcs, "-b")
    title("Recovered bcs")
    show()

    # figure(figsize=(20, 20))
    # subplot(2, 2, 1)
    # imshow(TH)
    # colorbar()
    # title("TH")
    # subplot(2, 2, 2)
    # imshow(inv(TH))
    # colorbar()
    # title("inv(TH)")
    # subplot(2, 2, 3)
    # imshow(K)
    # colorbar()
    # title("K")
    # subplot(2, 2, 4)
    # quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx], Ueff[1:2:end], Ueff[2:2:end])
    # title("Ueff")
    # show()

    
    #! Calculate and plot boundary displacements
    #! Not sure this will work because it didn't for Okada
    figure(figsize=(35,20))
    subplot(2, 4, 1)
    quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx], -Ueff[1:2:end], -Ueff[2:2:end])
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    title("-Ueff", fontsize=fontsize)

    subplot(2, 4, 5)
    plot(-Ueff[1:2:end], "-r", label="ux")
    plot(-Ueff[2:2:end], "-b", label="uy")
    gca().tick_params(labelsize=fontsize)
    legend()

    subplot(2, 4, 2)
    quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx], Uboundaryonly[1:2:end], Uboundaryonly[2:2:end])
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    title("T * Ueff", fontsize=fontsize)

    subplot(2, 4, 6)
    plot(Uboundaryonly[1:2:end], "-r", label="ux")
    plot(Uboundaryonly[2:2:end], "-b", label="uy")
    gca().tick_params(labelsize=fontsize)
    legend()

    subplot(2, 4, 3)
    quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx], -Ugravity[1:2:end], -Ugravity[2:2:end])
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    title("-Ugravity", fontsize=fontsize)

    subplot(2, 4, 7)
    plot(-Ugravity[1:2:end], "-r", label="ux")
    plot(-Ugravity[2:2:end], "-b", label="uy")
    gca().tick_params(labelsize=fontsize)
    legend()

    subplot(2, 4, 4)
    quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx], -Uall[1:2:end], -Uall[2:2:end])
    gca().set_aspect("equal")
    title("-Utotal = ", fontsize=fontsize)
    gca().tick_params(labelsize=fontsize)
    show()

    subplot(2, 4, 8)
    plot(-Uall[1:2:end], "-r", label="ux")
    plot(-Uall[2:2:end], "-b", label="uy")
    gca().tick_params(labelsize=fontsize)
    legend()


    @infiltrate
    return


    #! Interior displacements from boundaries (Ueff)
    # UinteriorBRTL, _ = constdispstress(slip2dispstress, x, y,
    #                                    els, BRTLidx,
    #                                    Ueff[1:2:end], Ueff[2:2:end],
    #                                    mu, nu)

    #! For each element centroid calculate the effect of a body force
    Uinteriorgravity = zeros(length(x), 2)
    Ugfield = zeros(length(x), 2)
    for i in 1:ntri
        xtri = mesh.point[1, mesh.cell[:, i]]
        ytri = mesh.point[2, mesh.cell[:, i]]
        triarea = trianglearea(xtri[1], ytri[1], xtri[2], ytri[3], xtri[3], ytri[3])
        tricentroid = [mean(xtri) mean(ytri)]
        Ugi, _ = kelvinUS(x, y, tricentroid[1], tricentroid[2], fx, fy, mu, nu)
        Uinteriorgravity += triarea .* Ugi
    end

    # figure()
    # quiver(x, y, Uinteriorgravity[:, 1], Uinteriorgravity[:, 2])
    # title("gravity")
    # show()

    # figure()
    # quiver(x, y, UinteriorBRTL[:, 1] .+ Uinteriorgravity[:, 1], UinteriorBRTL[:, 2] .+ Uinteriorgravity[:, 2])
    # title("boundaries + gravity")
    # show()

    figure(figsize=(30, 30))
    xmat = reshape(x, npts, npts)
    ymat = reshape(y, npts, npts)
    contourlevels = 20

    subplot(3, 3, 1)
    quiver(x, y, UinteriorBRTL[:, 1], UinteriorBRTL[:, 2])
    title("boundaries")

    subplot(3, 3, 2)
    mat = reshape(UinteriorBRTL[:, 1], npts, npts)
    contourf(xmat, ymat, mat, levels=contourlevels)
    colorbar()
    xlabel("x")
    ylabel("y")
    title("boundaries, ux")

    subplot(3, 3, 3)
    mat = reshape(UinteriorBRTL[:, 2], npts, npts)
    contourf(xmat, ymat, mat, levels=contourlevels)
    colorbar()
    xlabel("x")
    ylabel("y")
    title("boundaries, uy")

    subplot(3, 3, 4)
    quiver(x, y, Uinteriorgravity[:, 1],Uinteriorgravity[:, 2])
    title("gravity")

    subplot(3, 3, 5)
    mat = reshape(Uinteriorgravity[:, 1], npts, npts)
    contourf(xmat, ymat, mat, levels=contourlevels)
    colorbar()
    xlabel("x")
    ylabel("y")
    title("gravity, ux")

    subplot(3, 3, 6)
    mat = reshape(Uinteriorgravity[:, 2], npts, npts)
    contourf(xmat, ymat, mat, levels=contourlevels)
    colorbar()
    xlabel("x")
    ylabel("y")
    title("gravity, uy")


    subplot(3, 3, 7)
    quiver(x, y, UinteriorBRTL[:, 1] .- Uinteriorgravity[:, 1], UinteriorBRTL[:, 2] .- Uinteriorgravity[:, 2])
    title("boundaries + gravity")

    subplot(3, 3, 8)
    mat = reshape(UinteriorBRTL[:, 1] .- alpha .* Uinteriorgravity[:, 1], npts, npts)
    contourf(xmat, ymat, mat, levels=contourlevels)
    colorbar()
    xlabel("x")
    ylabel("y")
    title("boundaries + gravity, ux")

    subplot(3, 3, 9)
    mat = reshape(UinteriorBRTL[:, 2] .- alpha .* Uinteriorgravity[:, 2], npts, npts)
    contourf(xmat, ymat, mat, levels=contourlevels)
    colorbar()
    xlabel("x")
    ylabel("y")
    title("boundaries + gravity, uy")


    @infiltrate
    return



    # #! What I need is a full forward evaluation at the bottom coordinates only
    # Uboundarybottom, _ = constdispstress(slip2dispstress,
    #                                      els.xcenter[idx["bottom"]],
    #                                      els.ycenter[idx["bottom"]],
    #                                      els, collect(1:1:4*nels),
    #                                      Uboundarygravity[:, 1],
    #                                      Uboundarygravity[:, 2],
    #                                      mu, nu)

    # #! Forward evaluation of interior fields
    # Uvolumegravity, _ = constdispstress(slip2dispstress, x, y,
    #                                     els, collect(1:1:4*nels),
    #                                     Uboundarygravity[:, 1],
    #                                     Uboundarygravity[:, 2],
    #                                     mu, nu)
    # Uvolumebottom, _ = constdispstress(slip2dispstress, x, y,
    #                                    els, idx["bottom"],
    #                                    Uboundarybottom[:, 1], Uboundarybottom[:, 2],
    #                                    mu, nu)

    

    #! Plot meshes and fields
    figure(figsize=(20,20))
    nrows = 3
    ncols = 3
    fontsize = 20
    contour_levels = collect(LinRange(-200.0, 200.0, 100))
    contour_color = "white"
    contour_linewidth = 0.5
    markersize = 12

    #! Mesh
    subplot(nrows, ncols, 1)
    plotmesh(mesh)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("edge and area meshes", fontsize=fontsize)
    plot([x1, x2], [y1, y2], "-k", linewidth=2)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    #! Bottom boundary y-displacements
    subplot(nrows, ncols, 2)
    xmat = reshape(x, npts, npts)
    ymat = reshape(y, npts, npts)
    uymat = reshape(Uvolumegravity[:, 1], npts, npts)
    uymatdispgbottom = reshape(Uvolumebottom[:, 1], npts, npts)
    plot(xmat[:, 1], uymat[:, 1], "-r", label="Uvolumegravity", linewidth=2.0)
    plot(els.xcenter[idx["bottom"]], Uboundarybottom[:, 1], "-b", markersize=markersize, label="Uboundarybottom")
    plot(xmat[:, 1], 2 .* uymatdispgbottom[:, 1], "-k", markersize=markersize, label="Uvolumebottom")
    gca().tick_params(labelsize=fontsize)
    legend(fontsize=fontsize)
    title("x-displacements at bottom", fontsize=fontsize)

    #! Bottom boundary y-displacements
    subplot(nrows, ncols, 3)
    xmat = reshape(x, npts, npts)
    ymat = reshape(y, npts, npts)
    uymat = reshape(Uvolumegravity[:, 2], npts, npts)
    uymatdispgbottom = reshape(Uvolumebottom[:, 2], npts, npts)
    plot(xmat[:, 1], uymat[:, 1], "-r", label="Uvolumegravity", linewidth=2.0)
    plot(els.xcenter[idx["bottom"]], Uboundarybottom[:, 2], "-b", markersize=markersize, label="Uboundarybottom")
    plot(xmat[:, 1], 2 .* uymatdispgbottom[:, 1], "-k", markersize=markersize, label="Uvolumebottom")
    gca().tick_params(labelsize=fontsize)
    legend(fontsize=fontsize)
    title("y-displacements at bottom", fontsize=fontsize)

    #! Interior x-displacements
    subplot(nrows, ncols, 4)
    mat = Uvolumegravity[:, 1] #.- minimum(Uvolumegravity[:, 2]) 
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("ux (Uvolumegravity)", fontsize=fontsize)

    subplot(nrows, ncols, 5)
    mat = Uvolumebottom[:, 1] #.- minimum(Udispgbottom[:, 2])
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("ux (Uvolumebottom)", fontsize=fontsize)

    subplot(nrows, ncols, 6)
    mat1 = Uvolumegravity[:, 1] #.- minimum(Uvolumegravity[:, 2])
    mat2 = 2 .* Uvolumebottom[:, 1] #.- minimum(Uvolumebottom[:, 2])
    mat = mat1 .- mat2
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("ux (total)", fontsize=fontsize)

    #! Interior y-displacements
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
    title("uy (Uvolumegravity)", fontsize=fontsize)

    subplot(nrows, ncols, 8)
    mat = Uvolumebottom[:, 2] 
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("uy (Uvolumebottom)", fontsize=fontsize)

    subplot(nrows, ncols, 9)
    mat1 = Uvolumegravity[:, 2]
    mat2 = 2 .* Uvolumebottom[:, 2]
    mat = mat1 .- mat2
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    title("uy (total)", fontsize=fontsize)
    tight_layout()
    show()

    #! Quiver plot of total displacement field
    figure(figsize=(15, 15))
    xvec = Uvolumegravity[:, 1] .- 2 .* Uvolumebottom[:, 1] #.- minimum(Uvolumebottom[:, 2])
    yvec = Uvolumegravity[:, 2] .- 2 .* Uvolumebottom[:, 2] #.- minimum(Uvolumebottom[:, 2])
    quiver(x, y, xvec, yvec)
    show()
    # @infiltrate
end
gravitysquare()
