using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d


"""
    discretized_arc(thetastart, thetaend, radius, n_pts)

Generate regularly spaced eleemnts along an curved arc.
"""
function discretized_arc(thetastart, thetaend, radius, n_pts)
    # Create geometry of discretized arc
    thetarange = collect(LinRange(thetastart, thetaend, n_pts + 1))
    x = @. radius * cos(thetarange)
    y = @. radius * sin(thetarange)
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
    DEBUG = false
    muI = 1e10
    nuI = 0.25
    muII = 0.5 * muI
    nuII = 0.25
    p = muI / 1e3 # CS example
    nels = 600
    a = 0.5
    b = 1.0
    startangle = 180
    endangle = -startangle

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(startangle), deg2rad(endangle), a, nels)
    addelsez!(els, x1, y1, x2, y2, "aI")
    x1, y1, x2, y2 = discretized_arc(deg2rad(endangle), deg2rad(startangle), b, nels)
    addelsez!(els, x1, y1, x2, y2, "bI")
    # We want element K in bI to match exactly with element K in bII, so we
    # simply flip x1,y1 with x2,y2
    addelsez!(els, x2, y2, x1, y1, "bII")
    idx = getidxdict(els)

    if DEBUG
        figure()
        for n in ["aI", "bI"]
            plot(els.x1[idx[n]], els.y1[idx[n]])
            quiver(els.x1[idx[n]], els.y1[idx[n]], els.xnormal[idx[n]], els.ynormal[idx[n]])
        end
        gca().set_aspect("equal")
        title("region I boundaries and normals")

        figure()
        for n in ["bII"]
            plot(els.x1[idx[n]], els.y1[idx[n]])
            quiver(els.x1[idx[n]], els.y1[idx[n]], els.xnormal[idx[n]], els.ynormal[idx[n]])
        end
        gca().set_aspect("equal")
        title("region II boundaries and normals")
    end
    
    # Apply normal tractions everywhere and convert from radial to Cartesian
    xtracaI = zeros(length(idx["aI"]))
    ytracaI = zeros(length(idx["aI"]))
    for i in 1:length(idx["aI"]) # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtracaI[i], ytracaI[i] = els.rotmat[idx["aI"][i], :, :] * normalTractions
    end
    
    # Kernels and assembly
    TH = zeros(6*nels, 6*nels)

    # Region I materials
    T_aI_aI, H_aI_aI = PUTC(slip2dispstress, els, idx["aI"], idx["aI"], muI, nuI)    
    T_aI_bI, H_aI_bI = PUTC(slip2dispstress, els, idx["aI"], idx["bI"], muI, nuI)
    T_bI_aI, H_bI_aI = PUTC(slip2dispstress, els, idx["bI"], idx["aI"], muI, nuI)
    T_bI_bI, H_bI_bI = PUTC(slip2dispstress, els, idx["bI"], idx["bI"], muI, nuI)
    
    # Region II materials
    T_bII_bII, H_bII_bII = PUTC(slip2dispstress, els, idx["bII"], idx["bII"], muII, nuII)

    # Assemble BEM operator and boundary conditions
    # This row of equations enforces displacement equality at the boundary between regions.
    TH[1:(nels*2), 1:(nels*2)] = -T_bII_bII
    TH[1:(nels*2), (nels*2+1):(nels*4)] = T_bI_bI
    TH[1:(nels*2), (nels*4+1):(nels*6)] = T_bI_aI

    # This row of equations enforces traction continuity at the boundary between regions.
    TH[(nels*2+1):(nels*4), 1:(nels*2)] = H_bII_bII
    TH[(nels*2+1):(nels*4), (nels*2+1):(nels*4)] = H_bI_bI
    TH[(nels*2+1):(nels*4), (nels*4+1):(nels*6)] = H_bI_aI

    # This row of equations enforces the radial traction BC
    TH[(nels*4+1):(nels*6), (nels*2+1):(nels*4)] = H_aI_bI
    TH[(nels*4+1):(nels*6), (nels*4+1):(nels*6)] = H_aI_aI

    # Boundary conditions
    bcs = zeros(6*nels)
    bcs[(nels*4+1):(nels*6)] = interleave(xtracaI, ytracaI)

    # Simple diagonal preconditioner. Multiply every row by the inverse of its
    # diagonal entry so that the diagonal will be all ones.
    diag_entries = ones(6*nels)
    diag_entries[1:(6*nels)] = diag(TH)
    TH = TH ./ diag_entries
    bcs = bcs ./ diag_entries

    if DEBUG
        @show cond(TH)
        @show rank(TH)
        @show size(TH)
        figure()
        matshow(log10.(abs.(TH)))
        colorbar()
        figure()
        plot(els.x1[idx["aI"]], els.y1[idx["aI"]])
        quiver(els.x1[idx["aI"]], els.y1[idx["aI"]], bcs[(nels*4+1):2:end], bcs[(nels*4+2):2:end])
        title("r=a traction BCs")
    end
    
    # Solve BEM problem
    Ueff = TH \ bcs
    Ueffb2 = Ueff[1:1:(nels*2)]
    Ueffb1 = Ueff[(nels*2+1):1:(nels*4)]
    Ueffa1 = Ueff[(nels*4+1):1:(nels*6)]

    if DEBUG
        figure()
        plot(els.x1[idx["aI"]], els.y1[idx["aI"]])
        quiver(els.x1[idx["aI"]], els.y1[idx["aI"]], Ueffa1[1:2:end], Ueffa1[2:2:end])
        plot(els.x1[idx["bI"]], els.y1[idx["bI"]])
        quiver(els.x1[idx["bI"]], els.y1[idx["bI"]], Ueffb1[1:2:end], Ueffb1[2:2:end])
        title("region I Ueff")

        figure()
        plot(els.x1[idx["bII"]], els.y1[idx["bII"]])
        quiver(els.x1[idx["bII"]], els.y1[idx["bII"]], Ueffb2[1:2:end], Ueffb2[2:2:end])
        title("region II Ueff")

        # Effective displacements
        figure()
        plot(Ueff[1:2:end], ".r", label="ux")
        plot(Ueff[2:2:end], "+b", label="uy")
        legend()
        title("Ueff - whole vector")
    end
    
    # Forward line evaluation
    nprof = 30
    xprof = LinRange(0.51, 1.49, nprof)
    yprof = zeros(size(xprof))
    Uddm = zeros(nprof, 2)
    Sddm = zeros(nprof, 3)
    aidx = findall(x -> x <= b, xprof)
    bidx = findall(x -> x > b, xprof)
    UbII, SbII = constdispstress(slip2dispstress, xprof, yprof, els, idx["bII"],
                               Ueffb2[1:2:end], Ueffb2[2:2:end], muII, nuII)
    UbI, SbI = constdispstress(slip2dispstress, xprof, yprof, els, idx["bI"],
                               Ueffb1[1:2:end], Ueffb1[2:2:end], muI, nuI)
    UaI, SaI = constdispstress(slip2dispstress, xprof, yprof, els, idx["aI"],
                               Ueffa1[1:2:end], Ueffa1[2:2:end], muI, nuI)
    UI = UbI .+ UaI
    SI = SbI .+ SaI
    UII = UbII
    SII = SbII
    Uddm[aidx, :] = UI[aidx, :]
    Uddm[bidx, :] = UII[bidx, :]
    Sddm[aidx, :] = SI[aidx, :]
    Sddm[bidx, :] = SII[bidx, :]
    
    # Analytic solution
    nprofanalytic = 10000
    r = LinRange(0.5+1e-3, 1.50-1e-3, nprofanalytic)
    aidx = findall(x -> x <= b, r)
    bidx = findall(x -> x > b, r)
    pprime = (2*(1-nuI)*p*a^2/b^2) / (2*(1-nuI)+(muI/muII-1)*(1-a^2/b^2))
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
    text([0.96], [1.1], "I", horizontalalignment="center")
    text([1.04], [1.1], "II", horizontalalignment="center")
    plot(r, Stt ./ p, "-b", linewidth=linewidth, label="analytic")
    plot(xprof, Sddm[:, 2] ./ p, "ro", markersize=6, label="DDM")    
    xlabel(L"x \; / \; b")
    ylabel(L"\sigma_{yy} \; / \; p")
    xlim([0.5, 1.5])
    ylim([-0.3, 1.3])    
  
    legend()

    subplot(2, 1, 2)
    axvline(x=1, linewidth=1, linestyle=":", color="k")
    text([0.96], [0.1], "I", horizontalalignment="center")
    text([1.04], [0.1], "II", horizontalalignment="center")
    plot(r, Srr ./ p, "-b", linewidth=linewidth, label="analytic")
    plot(xprof, Sddm[:, 1] ./ p, "ro", markersize=6, label="DDM")    
    xlabel(L"x \; / \; b")
    ylabel(L"\sigma_{xx} \; / \; p")
    xlim([0.5, 1.5])
    ylim([-1.3, 0.3])    
    legend()
end
discmaterial()
