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
        fill(x, y, "lightgray", edgecolor="k", linewidth=0.25)
    end
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

            # TODO This is a really stupid area scaling.  Need to upgrade to Guassian 4-point.
            # Is there any good way to do this?
            _Ku *= triarea
            _Ks *= triarea

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
    gravitysquareparticular()

Experiments with gravity body force.
"""
function gravitysquareparticular()
    close("all")
    fontsize = 20
    mu = 3e10
    nu = 0.25
    rho = 2700
    g = 9.81
    nels = 20
    npts = 25
    L = 1e4
    # x, y = obsgrid(-L+1, -L+1, L-1, L-1, npts)
    x, y = obsgrid(-L+100, -L+100, L-100, L-100, npts)

    fx = 0.0
    fy = g * rho

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

    #! Kernels matrices (boundary -> boundary)
    Bidx = idx["B"]
    RTLidx = [idx["R"] ; idx["T"]; idx["L"]]
    BRTLidx = collect(1:1:els.endidx)
    T_pB_qBRTL, H_pB_qBRTL = PUTC(slip2dispstress, els, Bidx, BRTLidx, mu, nu)
    T_pRTL_qBRTL, H_pRTL_qBRTL = PUTC(slip2dispstress, els, RTLidx, BRTLidx, mu, nu)

    #! Particular solution technique with modified boundary conditoions
    # alpha = 7e-8
    # Tall = [T_pB_qBRTL ; T_pRTL_qBRTL]
    # bcs = zeros(2 * 4 * nels)
    # tempbc = @. alpha * -rho * g * (10000 - els.ycenter[idx["L"]])
    # bcs[42:2:80] = tempbc
    # tempbc = @. alpha * -rho * g * (10000 - els.ycenter[idx["R"]])
    # bcs[122:2:160] = tempbc
    # TH = [T_pB_qBRTL ; alpha .* H_pRTL_qBRTL]
    # Ueffparticular = inv(TH) * bcs
    # Uboundaryparticular = Tall * Ueffparticular
    

    # Modified Particular solution approach
    # Here I'm flipuding the left and right stresses and adding a top stress.
    # Just trying it out.
    alpha = 7e-8
    Tall = [T_pB_qBRTL ; T_pRTL_qBRTL]
    bcs = zeros(2 * 4 * nels)
    tempbc = @. alpha * -rho * g * (10000 - els.ycenter[idx["L"]])
    # bcs[122:2:160] = tempbc # Right hand side
    bcs[42:2:80] = tempbc # Left hand side
    tempbc = @. alpha * -rho * g * (10000 - els.ycenter[idx["R"]])
    bcs[122:2:160] = tempbc # Right hand side
    tempbc = @. alpha * -rho * g * 20000
    bcs[82:2:120] .= tempbc # Top


    # Try the Martel approach which has only traction BCs.  H matrix is horribly conditioned
    T_pBRTL_qBRTL, H_pBRTL_qBRTL = PUTC(slip2dispstress, els, BRTLidx, BRTLidx, mu, nu)


    figure()
    plot(bcs)
    show()
    @infiltrate
    return

    TH = [T_pB_qBRTL ; alpha .* H_pRTL_qBRTL]
    Ueffparticular = inv(TH) * bcs
    Uboundaryparticular = Tall * Ueffparticular



    #! Interior displacements from boundaries (Ueff)
    UinteriorBRTL, SinteriorBRTL = constdispstress(slip2dispstress, x, y, els, BRTLidx, Ueffparticular[1:2:end], Ueffparticular[2:2:end], mu, nu)
    figure()
    plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), UinteriorBRTL, SinteriorBRTL, "Particular integral method")
    
    #! Single quiver plot of interior displacements
    figure(figsize=(8, 8))
    quiver(x, y, UinteriorBRTL[:, 1], UinteriorBRTL[:, 2])
    xlabel("x (m)", fontsize=fontsize)
    ylabel("y (m)", fontsize=fontsize)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
    title("internal displacements", fontsize=fontsize)
end
gravitysquareparticular()
