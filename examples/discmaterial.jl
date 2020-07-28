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
    mu = 3e10
    nu = 0.25
    p = -1.0e5 # Applied radial pressure over arc
    nels = 360
    Ra = 0.5
    Rb = 1.0
    npts = 50
    x, y = obsgrid(-3, -3, 3, 3, npts)
    r = @. sqrt(x^2 + y^2)

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), Ra, nels)
    addelsez!(els, x1, y1, x2, y2, "a")
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), Rb, nels)
    addelsez!(els, x1, y1, x2, y2, "b")
    idx = getidxdict(els)

    # Apply normal tractions everywhere and convert from radial to Cartesian
    xtraca = zeros(length(idx["a"]))
    ytraca = zeros(length(idx["a"]))
    for i in 1:length(idx["a"]) # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtraca[i], ytraca[i] = els.rotmat[idx["a"][i], :, :] * normalTractions
    end
    
    # Kernels and assembly
    TH = zeros(6*nels, 6*nels)

    # Region 1 materials
    T_a1_a1, H_a1_a1 = PUTC(slip2dispstress, els, idx["a"], idx["a"], 2*mu, nu)    
    T_b1_a1, H_b1_a1 = PUTC(slip2dispstress, els, idx["b"], idx["a"], 2*mu, nu)
    T_b1_b1, H_b1_b1 = PUTC(slip2dispstress, els, idx["b"], idx["b"], 2*mu, nu)

    # Region 2 materials
    alpha = 1
    T_b2_b2, H_b2_b2 = PUTC(slip2dispstress, els, idx["b"], idx["b"], 1*mu, nu)
    TH[1:720, 1:720] = T_b2_b2
    TH[1:720, 721:1440] = T_b1_b1
    TH[1:720, 1441:2160] = T_b1_a1
    TH[721:1440, 1:720] = alpha .* H_b2_b2
    TH[721:1440, 721:1440] = alpha .* H_b1_b1
    TH[721:1440, 1441:2160] = alpha .* H_b1_a1
    TH[1441:2160, 1441:2160] = alpha .* H_a1_a1
    matshow(log10.(abs.(TH)))
    colorbar()
    bcs = zeros(6*nels)
    bcs[1441:2160] = interleave(xtraca, ytraca)
    @show cond(TH)
    
    # Solve BEM problem
    Ueff = TH \ bcs
    Ueffb2 = Ueff[1:1:720]
    Ueffb1 = Ueff[721:1:1440]
    Ueffa1 = Ueff[1441:1:2160]

    # Effective displacements
    # figure(figsize=(20, 10))
    # subplot(1, 2, 1)
    # plot(Ueff[1:2:end], ".r", label="ux")
    # plot(Ueff[2:2:end], "+b", label="uy")
    # legend()
    # title("Ueff - whole vector")

    # subplot(1, 2, 2)
    # plot(Ueffb2, ".r", label="ux, b2")
    # plot(Ueffb1, ".b", label="ux, b1")
    # plot(Ueffa1, ".g", label="ux, a1")
    # legend()
    # title("Ueff - subset selection")

    # figure(figsize=(20, 10))
    # subplot(1, 3, 1)
    # quiver(els.xcenter[idx["b"]], els.ycenter[idx["b"]],
    #        Ueffb2[1:2:end], Ueffb2[2:2:end])
    # gca().set_aspect("equal")
    # title("Ueff b2")

    # subplot(1, 3, 2)
    # quiver(els.xcenter[idx["b"]], els.ycenter[idx["b"]],
    #        Ueffb1[1:2:end], Ueffb1[2:2:end])
    # gca().set_aspect("equal")
    # title("Ueff b1")

    # subplot(1, 3, 3)
    # quiver(els.xcenter[idx["a"]], els.ycenter[idx["a"]],
    #        Ueffa1[1:2:end], Ueffa1[2:2:end])
    # gca().set_aspect("equal")
    # title("Ueff a1")

    # Forward line evaluation
    nprof = 100
    xprof = LinRange(0.51, 1.49, nprof)
    yprof = zeros(size(xprof))    
    Ub2, Sb2 = constdispstress(slip2dispstress, xprof, yprof, els,
                               idx["b"],
                               Ueffb2[1:2:end], Ueffb2[2:2:end], 1*mu, nu)
    Ub1, Sb1 = constdispstress(slip2dispstress, xprof, yprof, els,
                               idx["b"],
                               Ueffb1[1:2:end], Ueffb1[2:2:end], 2*mu, nu)
    Ua1, Sa1 = constdispstress(slip2dispstress, xprof, yprof, els,
                               idx["a"],
                               Ueffa1[1:2:end], Ueffa1[2:2:end], 2*mu, nu)

    # Analytic solution    
    figure(figsize=(6,6))
    subplot(2, 1, 1)
    plot(xprof, Sb2[:, 1] ./ p, "-g", label=L"\sigma_{yy}, b2")
    plot(xprof, Sb1[:, 1] ./ p, "-r", label=L"\sigma_{yy}, b1")
    plot(xprof, Sa1[:, 1] ./ p, "-b", label=L"\sigma_{yy}, a1")
    xlabel("x / b")
    ylabel(L"\sigma_{yy} / p")
    legend()

    subplot(2, 1, 2)
    plot(xprof, Sb2[:, 1] ./ p, "-g", label=L"\sigma_{xx}, b2")
    plot(xprof, Sb1[:, 1] ./ p, "-r", label=L"\sigma_{xx}, b1")
    plot(xprof, Sa1[:, 1] ./ p, "-b", label=L"\sigma_{xx}, a1")
    xlabel("x / b")
    ylabel(L"\sigma_{xx} / p")
    legend()
    
    # Forward volume evaluation
    # Ub2, Sb2 = constdispstress(slip2dispstress, x, y, els,
    #                            idx["b"],
    #                            Ueffb2[1:2:end], Ueffb2[2:2:end], 2*mu, nu)
    # Ub1, Sb1 = constdispstress(slip2dispstress, x, y, els,
    #                            idx["b"],
    #                            Ueffb1[1:2:end], Ueffb1[2:2:end], mu, nu)
    # Ua1, Sa1 = constdispstress(slip2dispstress, x, y, els,
    #                            idx["a"],
    #                            Ueffa1[1:2:end], Ueffa1[2:2:end], mu, nu)
        
    # # Summary figure
    # figure(figsize=(20,10))
    # nrows = 3
    # ncols = 4
    
    # # Volume solutions
    # circle_subplot(nrows, ncols, 1, els, x, y, Ub2[:, 1], npts, "ux (b2)")
    # circle_subplot(nrows, ncols, 2, els, x, y, Ub2[:, 2], npts, "uy (b2)")
    # circle_subplot(nrows, ncols, 5, els, x, y, Ub1[:, 1], npts, "ux (b1)")
    # circle_subplot(nrows, ncols, 6, els, x, y, Ub1[:, 2], npts, "uy (b1)")
    # circle_subplot(nrows, ncols, 9, els, x, y, Ua1[:, 1], npts, "ux (a1)")
    # circle_subplot(nrows, ncols, 10, els, x, y, Ua1[:, 2], npts, "uy (a1)")
end
discmaterial()
