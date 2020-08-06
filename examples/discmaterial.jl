using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
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
function circle_subplot(nrows, ncols, plotidx, els, x, y, mat, npts, title_string)
    contour_levels = 100
    contour_color = "white"
    contour_linewidth = 0.5

    subplot(nrows, ncols, plotidx)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.05, pad=0.05, extend = "both")
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    xlabel("x (m)")
    ylabel("y (m)")
    title(title_string)
    plotelements(els)
    gca().set_aspect("equal")
end


"""
    discmaterial()

An attmept at the Crouch and Starfield level annulus solution
with varying material properties
"""
function discmaterial()
    close("all")
    mu1 = 1e10
    nu1 = 0.25
    mu2 = 0.5 * mu1
    nu2 = 0.25
    p = mu1 / 1e3 # CS example
    nels = 600
    a = 0.5
    b = 1.0
    npts = 50
    x, y = obsgrid(-3, -3, 3, 3, npts)
    r = @. sqrt(x^2 + y^2)

    start_angle = 180
    end_angle = -start_angle

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(start_angle), deg2rad(end_angle), a, nels)
    addelsez!(els, x1, y1, x2, y2, "a")
    x1, y1, x2, y2 = discretized_arc(deg2rad(end_angle), deg2rad(start_angle), b, nels)
    addelsez!(els, x1, y1, x2, y2, "b_I")
    # We want element K in b_I to match exactly with element K in b_II, so we
    # simply flip x1,y1 with x2,y2
    addelsez!(els, x2, y2, x1, y1, "b_II")
    idx = getidxdict(els)

    figure()
    for n in ["a", "b_I"]
        plot(els.x1[idx[n]], els.y1[idx[n]])
        quiver(els.x1[idx[n]], els.y1[idx[n]], els.xnormal[idx[n]], els.ynormal[idx[n]])
    end
    title("region I boundaries and normals")

    figure()
    for n in ["b_II"]
        plot(els.x1[idx[n]], els.y1[idx[n]])
        quiver(els.x1[idx[n]], els.y1[idx[n]], els.xnormal[idx[n]], els.ynormal[idx[n]])
    end
    title("region II boundaries and normals")

    # Apply normal tractions everywhere and convert from radial to Cartesian
    xtraca = zeros(length(idx["a"]))
    ytraca = zeros(length(idx["a"]))
    for i in 1:length(idx["a"]) # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtraca[i], ytraca[i] = els.rotmat[idx["a"][i], :, :] * normalTractions
    end
    
    # Kernels and assembly
    # TH = zeros(6*nels+4, 6*nels)
    TH = zeros(6*nels, 6*nels)

    # Region 1 materials
    T_a1_a1, H_a1_a1 = PUTC(slip2dispstress, els, idx["a"], idx["a"], mu1, nu1)    
    T_a1_b1, H_a1_b1 = PUTC(slip2dispstress, els, idx["a"], idx["b_I"], mu1, nu1)
    T_b1_a1, H_b1_a1 = PUTC(slip2dispstress, els, idx["b_I"], idx["a"], mu1, nu1)
    T_b1_b1, H_b1_b1 = PUTC(slip2dispstress, els, idx["b_I"], idx["b_I"], mu1, nu1)
    
    # Region 2 materials
    T_b2_b2, H_b2_b2 = PUTC(slip2dispstress, els, idx["b_II"], idx["b_II"], mu2, nu2)

    # Assemble BEM operator and boundary conditions
    # This row of equations enforces displacement equality at the boundary between regions.
    TH[1:(nels*2), 1:(nels*2)] = -T_b2_b2
    TH[1:(nels*2), (nels*2+1):(nels*4)] = T_b1_b1
    TH[1:(nels*2), (nels*4+1):(nels*6)] = T_b1_a1

    # This row of equations enforces traction continuity at the boundary between regions.
    TH[(nels*2+1):(nels*4), 1:(nels*2)] = H_b2_b2
    TH[(nels*2+1):(nels*4), (nels*2+1):(nels*4)] = H_b1_b1
    TH[(nels*2+1):(nels*4), (nels*4+1):(nels*6)] = H_b1_a1

    # This row of equations enforces the radial traction BC
    TH[(nels*4+1):(nels*6), (nels*2+1):(nels*4)] = H_a1_b1
    TH[(nels*4+1):(nels*6), (nels*4+1):(nels*6)] = H_a1_a1

    # Avoid rigid body translations and rotations. TODO: I don't think this is
    # quite right. Possibly unnecessary
    # TH[nels*6+1, nels*2+1] = 1.0
    # TH[nels*6+2, nels*2+2] = 1.0
    # TH[nels*6+3, nels*4+3] = 1.0
    # TH[nels*6+4, nels*4+4] = 1.0

    bcs = zeros(6*nels)
    bcs[(nels*4+1):(nels*6)] = interleave(xtraca, ytraca)

    # Simple diagonal preconditioner. Multiply every row by the inverse of its
    # diagonal entry so that the diagonal will be all ones.
    diag_entries = ones(6*nels)
    diag_entries[1:(6*nels)] = diag(TH)
    TH = TH ./ diag_entries
    bcs = bcs ./ diag_entries

    @show cond(TH)
    @show rank(TH)
    @show size(TH)
    
    figure()
    matshow(log10.(abs.(TH)))
    colorbar()

    figure()
    plot(els.x1[idx["a"]], els.y1[idx["a"]])
    quiver(els.x1[idx["a"]], els.y1[idx["a"]], bcs[(nels*4+1):2:end], bcs[(nels*4+2):2:end])
    title("r=a traction BCs")
    
    # Solve BEM problem
    Ueff = TH \ bcs # Looks more reasonable but is it?
    # Ueff = inv(TH) * bcs # Horrible solution

    Ueffb2 = Ueff[1:1:(nels*2)]
    Ueffb1 = Ueff[(nels*2+1):1:(nels*4)]
    Ueffa1 = Ueff[(nels*4+1):1:(nels*6)]

    figure()
    plot(els.x1[idx["a"]], els.y1[idx["a"]])
    quiver(els.x1[idx["a"]], els.y1[idx["a"]], Ueffa1[1:2:end], Ueffa1[2:2:end])
    plot(els.x1[idx["b_I"]], els.y1[idx["b_I"]])
    quiver(els.x1[idx["b_I"]], els.y1[idx["b_I"]], Ueffb1[1:2:end], Ueffb1[2:2:end])
    title("region I Ueff")

    figure()
    plot(els.x1[idx["b_II"]], els.y1[idx["b_II"]])
    quiver(els.x1[idx["b_II"]], els.y1[idx["b_II"]], Ueffb2[1:2:end], Ueffb2[2:2:end])
    title("region II Ueff")

    # Effective displacements
    figure()
    plot(Ueff[1:2:end], ".r", label="ux")
    plot(Ueff[2:2:end], "+b", label="uy")
    legend()
    title("Ueff - whole vector")
    
    # Forward line evaluation
    nprof = 30
    xprof = LinRange(0.51, 1.49, nprof)
    yprof = zeros(size(xprof))
    Uddm = zeros(nprof, 2)
    Sddm = zeros(nprof, 3)
    aidx = findall(x -> x <= b, xprof)
    bidx = findall(x -> x > b, xprof)
    Ub2, Sb2 = constdispstress(slip2dispstress, xprof, yprof, els, idx["b_II"],
                               Ueffb2[1:2:end], Ueffb2[2:2:end], mu2, nu2)
    Ub1, Sb1 = constdispstress(slip2dispstress, xprof, yprof, els, idx["b_I"],
                               Ueffb1[1:2:end], Ueffb1[2:2:end], mu1, nu1)
    Ua1, Sa1 = constdispstress(slip2dispstress, xprof, yprof, els, idx["a"],
                               Ueffa1[1:2:end], Ueffa1[2:2:end], mu1, nu1)
    U1 = Ub1 .+ Ua1
    S1 = Sb1 .+ Sa1
    U2 = Ub2
    S2 = Sb2
    Uddm[aidx, :] = U1[aidx, :]
    Uddm[bidx, :] = U2[bidx, :]
    Sddm[aidx, :] = S1[aidx, :]
    Sddm[bidx, :] = S2[bidx, :]
    

    # Analytic solution
    nprofanalytic = 10000
    r = LinRange(0.5+1e-3, 1.50-1e-3, nprofanalytic)
    aidx = findall(x -> x <= b, r)
    bidx = findall(x -> x > b, r)
    pprime = (2*(1-nu1)*p*a^2/b^2) / (2*(1-nu1)+(mu1/mu2-1)*(1-a^2/b^2))
    a2b2 = a^2/b^2
    Srr1 = @. 1/(1-a2b2) * (p*a2b2-pprime - (p-pprime)*a^2/r^2)
    Stt1 = @. 1/(1-a2b2) * (p*a2b2-pprime + (p-pprime)*a^2/r^2)
    Srr2 = @. -pprime * b^2 / r^2
    Stt2 = @. pprime * b^2 / r^2
    Srr = zeros(size(Srr1))
    Stt = zeros(size(Stt1))
    Srr[aidx] = Srr1[aidx]
    Srr[bidx] = Srr2[bidx]
    Stt[aidx] = Stt1[aidx]
    Stt[bidx] = Stt2[bidx]
    
    # Graphically compare analytic and BEM solutions
    linewidth = 2.0
    figure(figsize=(7, 7))
    subplot(2, 1, 1)
    axvline(x=1, linewidth=1, linestyle=":", color="k")
    plot(r, Stt ./ p, "-b", linewidth=linewidth, label="analytic")
    plot(xprof, Sddm[:, 2] ./ p, "ro", markersize=6, label="DDM")    
    xlabel(L"x \; / \; b")
    ylabel(L"\sigma_{yy} \; / \; p")
    xlim([0.5, 1.5])
    legend()

    subplot(2, 1, 2)
    axvline(x=1, linewidth=1, linestyle=":", color="k")
    plot(r, Srr ./ p, "-b", linewidth=linewidth, label="analytic")
    plot(xprof, Sddm[:, 1] ./ p, "ro", markersize=6, label="DDM")    
    xlabel(L"x \; / \; b")
    ylabel(L"\sigma_{xx} \; / \; p")
    xlim([0.5, 1.5])
    legend()
    # savefig("yes.pdf")
end
discmaterial()
