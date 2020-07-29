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
    discmaterialcompress()

An attmept compressing a disc with an internal circle of higher modulus
"""
function discmaterialcompress()
    close("all")
    mu1 = 1e10
    nu1 = 0.25
    mu2 = 1.0 * mu1
    nu2 = 0.25
    p = mu1 / 1e3 # CS example
    du = 1e-5
    nels = 360
    a = 0.5
    b = 1.0
    npts = 50
    x, y = obsgrid(-1, -1, 1, 1, npts)
    r = @. sqrt(x^2 + y^2)

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), a, nels)
    addelsez!(els, x1, y1, x2, y2, "a")
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), b, nels)
    addelsez!(els, x1, y1, x2, y2, "b")
    idx = getidxdict(els)

    # Apply normal tractions everywhere and convert from radial to Cartesian
    xtracb = zeros(length(idx["b"]))
    ytracb = zeros(length(idx["b"]))
    for i in 1:length(idx["a"]) # Calcuate the x and y components of the tractions
        normaldisp = [0; du] # Pressure in fault normal component only.
        xtracb[i], ytracb[i] = els.rotmat[idx["b"][i], :, :] * normaldisp
    end
    
    # Kernels and assembly
    TH = zeros(6*nels, 6*nels)

    # Region 1 materials
    T_a1_a1, H_a1_a1 = PUTC(slip2dispstress, els, idx["a"], idx["a"], mu1, nu1)    
    T_a1_b1, H_a1_b1 = PUTC(slip2dispstress, els, idx["a"], idx["b"], mu1, nu1)
    T_b1_a1, H_b1_a1 = PUTC(slip2dispstress, els, idx["b"], idx["a"], mu1, nu1)
    T_b1_b1, H_b1_b1 = PUTC(slip2dispstress, els, idx["b"], idx["b"], mu1, nu1)
    
    # Region 2 materials
    T_a2_a2, H_a2_a2 = PUTC(slip2dispstress, els, idx["a"], idx["a"], mu1, nu1)

    # Assemble BEM operator and boundary conditions
    alpha = 1e0
    TH[1:720, 1:720] = T_b1_b1
    TH[1:720, 721:1440] = T_b1_a1
    TH[1:720, 1441:2160] .= 0
    TH[721:1440, 1:720] = T_a1_b1
    TH[721:1440, 721:1440] = T_a1_a1
    TH[721:1440, 1441:2160] = -T_a2_a2
    TH[1441:2160, 1:720] = H_a1_b1
    TH[1441:2160, 721:1440] = H_a1_a1
    TH[1441:2160, 1441:2160] = H_a2_a2

    @show cond(TH)
    @show rank(TH)
    bcs = zeros(6*nels)
    bcs[1:720] = interleave(xtracb, ytracb)
    matshow(log10.(abs.(TH)))
    colorbar()
    
    # Solve BEM problem
    Ueff = TH \ bcs # Looks more reasonable but is it?
    Ueffb1 = Ueff[1:1:720]
    Ueffa1 = Ueff[721:1:1440]
    Ueffa2 = Ueff[1441:1:2160]

    Ueffb = T_b1_b1 \ interleave(xtracb, ytracb)
    @show cond(T_b1_b1)
    @show rank(T_b1_b1)
    
    # Effective displacements
    figure()
    plot(Ueff[1:2:end], ".r", label="ux")
    plot(Ueff[2:2:end], "+b", label="uy")
    legend()
    title("Ueff - whole vector")
    
    # Forward line evaluation
    nprof = 70
    xprof = LinRange(0.51, 1.49, nprof)
    yprof = zeros(size(xprof))
    Ub1, Sb1 = constdispstress(slip2dispstress, x, y, els, idx["b"],
                               Ueffb1[1:2:end], Ueffb1[2:2:end], mu2, nu2)
    Ua1, Sa1 = constdispstress(slip2dispstress, x, y, els, idx["a"],
                               Ueffa1[1:2:end], Ueffa1[2:2:end], mu1, nu1)
    Ua2, Sa2 = constdispstress(slip2dispstress, x, y, els, idx["a"],
                               Ueffa2[1:2:end], Ueffa2[2:2:end], mu1, nu1)
    Ub, Sb = constdispstress(slip2dispstress, x, y, els, idx["b"],
                             Ueffb[1:2:end], Ueffb[2:2:end], mu1, nu1)

    
    nrows = 4
    ncols = 3
    figure(figsize=(10, 10))
    circle_subplot(nrows, ncols, 1, els, x, y, Ub1[:, 1], npts, "bI")
    circle_subplot(nrows, ncols, 2, els, x, y, Ub1[:, 2], npts, "bI")
    circle_subplot(nrows, ncols, 3, els, x, y, sqrt.(Ub1[:, 1].^2 + Ub1[:, 2].^2), npts, "bI")

    circle_subplot(nrows, ncols, 4, els, x, y, Ua1[:, 1], npts, "aI")
    circle_subplot(nrows, ncols, 5, els, x, y, Ua1[:, 2], npts, "aI")
    circle_subplot(nrows, ncols, 6, els, x, y, sqrt.(Ua1[:, 1].^2 + Ua1[:, 2].^2), npts, "aI")

    circle_subplot(nrows, ncols, 7, els, x, y, Ua2[:, 1], npts, "aII")
    circle_subplot(nrows, ncols, 8, els, x, y, Ua2[:, 2], npts, "aII")
    circle_subplot(nrows, ncols, 9, els, x, y, sqrt.(Ua2[:, 1].^2 + Ua2[:, 2].^2), npts, "aII")    

    circle_subplot(nrows, ncols, 10, els, x, y, Ub[:, 1], npts, "b all alone")
    circle_subplot(nrows, ncols, 11, els, x, y, Ub[:, 2], npts, "b all alone")
    circle_subplot(nrows, ncols, 12, els, x, y, sqrt.(Ub[:, 1].^2 + Ub[:, 2].^2), npts, "b all alone")    

    
end
discmaterialcompress()
