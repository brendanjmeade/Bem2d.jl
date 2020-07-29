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
    nels = 360
    a = 0.5
    b = 1.0
    npts = 50
    x, y = obsgrid(-3, -3, 3, 3, npts)
    r = @. sqrt(x^2 + y^2)

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), a, nels)
    addelsez!(els, x1, y1, x2, y2, "a")
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), b, nels)
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
    T_a1_a1, H_a1_a1 = PUTC(slip2dispstress, els, idx["a"], idx["a"], mu1, nu1)    
    T_a1_b1, H_a1_b1 = PUTC(slip2dispstress, els, idx["a"], idx["b"], mu1, nu1)
    T_b1_a1, H_b1_a1 = PUTC(slip2dispstress, els, idx["b"], idx["a"], mu1, nu1)
    T_b1_b1, H_b1_b1 = PUTC(slip2dispstress, els, idx["b"], idx["b"], mu1, nu1)

    # Region 2 materials
    T_b2_b2, H_b2_b2 = PUTC(slip2dispstress, els, idx["b"], idx["b"], mu2, nu2)

    # Assemble BEM operator and boundary conditions
    alpha = 1e-10
    TH[1:720, 1:720] = T_b2_b2
    TH[1:720, 721:1440] = T_b1_b1
    TH[1:720, 1441:2160] = T_b1_a1
    TH[721:1440, 1:720] = alpha .* H_b2_b2
    TH[721:1440, 721:1440] = alpha .* H_b1_b1
    TH[721:1440, 1441:2160] = alpha .* H_b1_a1
    TH[1441:2160, 721:1440] = alpha .* H_a1_b1 # Was I missing this before?
    TH[1441:2160, 1441:2160] = alpha .* H_a1_a1
    @show cond(TH)
    @show rank(TH)
    bcs = zeros(6*nels)
    bcs[1441:2160] = interleave(xtraca, ytraca)
    matshow(log10.(abs.(TH)))
    colorbar()
    
    # Solve BEM problem
    Ueff = TH \ bcs # Looks more reasonable but is it?
    # Ueff = inv(TH) * bcs # Horrible solution

    Ueffb2 = Ueff[1:1:720]
    Ueffb1 = Ueff[721:1:1440]
    Ueffa1 = Ueff[1441:1:2160]

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
    # aidx = findall(x -> x <= b, xprof)
    # bidx = findall(x -> x > b, xprof)
    Ub2, Sb2 = constdispstress(slip2dispstress, xprof, yprof, els, idx["b"],
                               Ueffb2[1:2:end], Ueffb2[2:2:end], mu2, nu2)
    Ub1, Sb1 = constdispstress(slip2dispstress, xprof, yprof, els, idx["b"],
                               Ueffb1[1:2:end], Ueffb1[2:2:end], mu1, nu1)
    Ua1, Sa1 = constdispstress(slip2dispstress, xprof, yprof, els, idx["a"],
                               Ueffa1[1:2:end], Ueffa1[2:2:end], mu1, nu1)

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
    linewidth = 1.0
    figure(figsize=(8, 8))
    subplot(2, 1, 1)
    plot(r, Stt ./ p, "-k", linewidth=linewidth, label="analytic")
    plot(xprof, Sb2[:, 1] ./ p, ".g", label=L"\sigma_{yy}, b^{II}")
    plot(xprof, Sb1[:, 1] ./ p, ".r", label=L"\sigma_{yy}, b^{I}")
    plot(xprof, Sa1[:, 1] ./ p, ".c", label=L"\sigma_{yy}, a^{I}")
    xlabel(L"x \; / \; b")
    ylabel(L"\sigma_{yy} \; / \; p")
    xlim([0.5, 1.5])
    # ylim([-1.2, 1.2])
    legend()

    subplot(2, 1, 2)
    plot(r, Srr ./ p, "-k", linewidth=linewidth, label="analytic")
    plot(xprof, Sb2[:, 1] ./ p, ".g", label=L"\sigma_{xx}, b^{II}")
    plot(xprof, Sb1[:, 1] ./ p, ".r", label=L"\sigma_{xx}, b^{I}")
    plot(xprof, Sa1[:, 1] ./ p, ".c", label=L"\sigma_{xx}, a^{I}")
    xlabel(L"x \; / \; b")
    ylabel(L"\sigma_{xx} \; / \; p")
    xlim([0.5, 1.5])
    # ylim([-1.2, 1.2])
    legend()
end
discmaterial()
