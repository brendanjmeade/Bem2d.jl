using Revise
using PyPlot
using PyCall
using Statistics
using LinearAlgebra
using Infiltrator
using Bem2d


"""
    plotbcsUeff(els, titlestring)

Plot boundary conditions and Ueff
"""
function plotbcsUeff(els, bcs, Ueff, suptitlestring)
    figure(figsize=(20, 20))
    # subplot(2, 2, 1)
    # quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx],
    #        bcs[1:2:end], bcs[2:2:end])
    # gca().set_aspect("equal")
    # title("bcs")
    
    subplot(2, 2, 1)
    plot(bcs[1:2:end], "rx", label="x")
    plot(bcs[2:2:end], "b+", label="y")
    legend()
    title("bcs")

    subplot(2, 2, 3)
    quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx],
           Ueff[1:2:end], Ueff[2:2:end])
    gca().set_aspect("equal")
    title("Ueff")

    subplot(2, 2, 4)
    plot(Ueff[1:2:end], "rx", label="x")
    plot(Ueff[2:2:end], "b+", label="y")
    legend()
    title("Ueff")

    suptitle(suptitlestring)
end        


"""
    plotgeometry(els, titlestring)

els structure geometry plotting for diagnostics
"""
function plotgeometry(els, titlestring)
    scale = 500
    figure(figsize=(15,10))
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-g")
        quiver(els.xcenter[i], els.ycenter[i],
               els.xnormal[i], els.ynormal[i],
               width=1e-3, scale=1e2)
        text(els.xcenter[i] + scale * els.xnormal[i],
             els.ycenter[i] + scale * els.ynormal[i],
             string(i), fontsize=10)
    end
    title(titlestring)
    gca().set_aspect("equal")
end        


"""
    flipud()

Just an alias to the 1-linear that replicates matlab's matrix flipud
Taken from:
https://cheatsheets.quantecon.org/#manipulating-vectors-and-matrices
"""
function flipud(mat)
    mat = reverse(mat, dims = 1)
end


"""
    layeredboxfault()

Comparing half-space and dislocation in a box solutions
"""
function layeredboxfault()
    close("all")
    PLOTGEOMETRY = false
    mu = 3e10
    nu = 0.25
    mu1 = 3e10
    nu1 = 0.25
    mu2 = 3e10
    nu2 = 0.25
    npts = 50
    offset = 100 # meters
    
    # Element geometries and data structures for the box case
    elsbox = Elements(Int(1e5))
    nfault = 1
    boxB = -30e3
    boxR = 30e3
    boxT = 0
    boxL = -30e3    
    nside = 20
    x1, y1, x2, y2 = discretizedline(boxL, boxB, boxR, boxB, nside) # Bottom
    addelsez!(elsbox, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(boxR, boxB, boxR, boxT, nside) # Right hand side
    addelsez!(elsbox, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(boxR, boxT, boxL, boxT, nside) # Top
    addelsez!(elsbox, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(boxL, boxT, boxL, boxB, nside) # Left hand side
    addelsez!(elsbox, x1, y1, x2, y2, "L")
    idxbox = getidxdict(elsbox)
    PLOTGEOMETRY && plotgeometry(elsbox, "Box boundaries and normals")

    # Box BEM problem
    T_B_BRTL, _ = PUTC(slip2dispstress, elsbox, idxbox["B"], 1:elsbox.endidx, mu, nu)
    _, H_RTL_BRTL = PUTC(slip2dispstress, elsbox, [idxbox["R"]; idxbox["T"]; idxbox["L"]],
                         1:elsbox.endidx, mu, nu)
    bcsbox = zeros(2*elsbox.endidx)
    bcsbox[2*51] = 1.0
    bcsbox[2*50] = 1.0    
    Ueffbox = [T_B_BRTL ; H_RTL_BRTL] \ bcsbox
    # plotbcsUeff(elsbox, bcsbox, Ueffbox, "Box")
    xgrid, ygrid = obsgrid(boxL+offset, boxB+offset, boxR-offset, boxT-offset, npts)
    U1, S1 = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, 1:elsbox.endidx,
                             Ueffbox[1:2:end], Ueffbox[2:2:end], mu, nu)
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               U1, S1, "Homogenous box")
    
    # Element geometries and data structures for the layered box case
    elslayer = Elements(Int(1e5))
    l1B = -15e3
    l1R = 30e3
    l1T = 0
    l1L = -30e3
    l2B = -30e3
    l2R = 30e3
    l2T = -15e3
    l2L = -30e3
    nside = 20
    x1, y1, x2, y2 = discretizedline(l1L, l1B, l1R, l1B, nside) # Layer 1 bottom
    addelsez!(elslayer, x1, y1, x2, y2, "B1")
    x1, y1, x2, y2 = discretizedline(l1R, l1B, l1R, l1T, nside) # Layer 1 right
    addelsez!(elslayer, x1, y1, x2, y2, "R1")
    x1, y1, x2, y2 = discretizedline(l1R, l1T, l1L, l1T, nside) # Layer 1 top
    addelsez!(elslayer, x1, y1, x2, y2, "T1")
    x1, y1, x2, y2 = discretizedline(l1L, l1T, l1L, l1B, nside) # Layer 1 left
    addelsez!(elslayer, x1, y1, x2, y2, "L1")
    x1, y1, x2, y2 = discretizedline(l2L, l2B, l2R, l2B, nside) # Layer 2 bottom
    addelsez!(elslayer, x1, y1, x2, y2, "B2")
    x1, y1, x2, y2 = discretizedline(l2R, l2B, l2R, l2T, nside) # Layer 2 right
    addelsez!(elslayer, x1, y1, x2, y2, "R2")

    # Clockwise ordering
    x1, y1, x2, y2 = discretizedline(l2R, l2T, l2L, l2T, nside) # Layer 2 top
    addelsez!(elslayer, x1, y1, x2, y2, "T2")

    # Ben's reording trick
    # x1, y1, x2, y2 = discretizedline(l2R, l2T, l2L, l2T, nside) # Layer 2 top
    # addelsez!(elslayer, x2, y2, x1, y1, "T2")

    x1, y1, x2, y2 = discretizedline(l2L, l2T, l2L, l2B, nside) # Layer 2 left
    addelsez!(elslayer, x1, y1, x2, y2, "L2")
    idxlayer = getidxdict(elslayer)
    plotgeometry(elslayer, "Layer boundaries and normals")

    # Build design matrix
    C1idx = [idxlayer["B1"] ; idxlayer["R1"] ; idxlayer["T1"] ; idxlayer["L1"]]
    C2idx = [idxlayer["B2"] ; idxlayer["R2"] ; idxlayer["T2"] ; idxlayer["L2"]]
    
    # Boundaries shared by C1 and C2
    T_B1_C1, H_B1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["B1"], C1idx, mu, nu)
    T_T2_C2, H_T2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["T2"], C2idx, mu, nu)

    # C1 only boundaries
    T_R1_C1, H_R1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["R1"], C1idx, mu, nu)
    T_T1_C1, H_T1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["T1"], C1idx, mu, nu)
    T_L1_C1, H_L1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["L1"], C1idx, mu, nu)

    # C2 only boundaries
    T_B2_C2, H_B2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["B2"], C2idx, mu, nu)
    T_R2_C2, H_R2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["R2"], C2idx, mu, nu)
    T_L2_C2, H_L2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["L2"], C2idx, mu, nu)

    # Place submatrices into larger matrix
    TH = zeros(2*elslayer.endidx, 2*elslayer.endidx)
    TH[1:2*nside, 1:8*nside] = T_B1_C1 # Displacements at B1 due to C1
    TH[1:2*nside, 8*nside+1:16*nside] = -T_T2_C2 # Displacements at T2 due to C2
    TH[2*nside+1:4*nside, 1:8*nside] = H_B1_C1 # Tractions at B1 due to C1
    TH[2*nside+1:4*nside, 8*nside+1:16*nside] = H_T2_C2 # Tractions at T2 due to C2
    TH[4*nside+1:6*nside, 1:8*nside] = H_R1_C1 # Tractions at R1 due to C1
    TH[6*nside+1:8*nside, 1:8*nside] = H_T1_C1 # Tractions at T1 due to C1
    TH[8*nside+1:10*nside, 1:8*nside] = H_L1_C1 # Tractions at L1 due to C1
    TH[10*nside+1:12*nside, 8*nside+1:16*nside] = T_B2_C2 # Displacments at B2 due to C2
    TH[12*nside+1:14*nside, 8*nside+1:16*nside] = H_R2_C2 # Tractions at R2 due to C2
    TH[14*nside+1:16*nside, 8*nside+1:16*nside] = H_L2_C2 # Tractions at L2 due to C2

    bcslayer = zeros(2*elslayer.endidx)
    bcslayer[140] = 1.0
    bcslayer[142] = 1.0
    Uefflayer = TH \ bcslayer
    plotbcsUeff(elslayer, bcslayer, Uefflayer, "Two-layer")
    UeffT2 = Uefflayer[1:elslayer.endidx]
    UeffB2 = Uefflayer[elslayer.endidx+1:end]
            
    # Volume visualization for two domain case
    Tidx = findall(x -> x > -15e3, ygrid)
    Bidx = findall(x -> x < -15e3, ygrid)
    Tx = xgrid[Tidx]
    Ty = ygrid[Tidx]
    Bx = xgrid[Bidx]
    By = ygrid[Bidx]

    UT2, ST2 = constdispstress(slip2dispstress, Tx, Ty, elslayer,
                               C1idx,
                               UeffT2[1:2:end], UeffT2[2:2:end], mu1, nu1)
    UB2, SB2 = constdispstress(slip2dispstress, Bx, By, elslayer,
                               C2idx,
                               UeffB2[1:2:end], UeffB2[2:2:end], mu2, nu2)
    # UT2F, ST2F = constdispstress(slip2dispstress, Tx, Ty, els2, idx2["F"],
    #                              Fslip[1:2:end], Fslip[2:2:end], mu1, nu1)
    # UB2F, SB2F = constdispstress(slip2dispstress, Bx, By, els2, idx2["F"],
    #                              Fslip[1:2:end], Fslip[2:2:end], mu2, nu2)    
    # UT2 = UT2 .+ UT2F
    # UB2 = UB2 .+ UB2F
    # ST2 = ST2 .+ ST2F
    # SB2 = SB2 .+ SB2F
    U2 = [UB2 ; UT2] # Note B, T ordering
    S2 = [SB2 ; ST2] # Note B, T ordering
    plotfields(elslayer, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               U2, S2, "two-layer")
    # plotfields(els2, reshape(xgrid1, npts, npts), reshape(ygrid1, npts, npts),
    #            U2 .- U1, S2 .- S1, "residuals")    

end
layeredboxfault()
