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
    subplot(2, 2, 1)
    quiver(els.xcenter[1:els.endidx], els.ycenter[1:els.endidx],
           bcs[1:2:end], bcs[2:2:end])
    gca().set_aspect("equal")
    title("bcs")
    
    subplot(2, 2, 2)
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
             string(i), fontsize=5)
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
    dislocationlayers()

Comparing half-space and dislocaiton in a box solutions
"""
function dislocationlayers()
    close("all")
    PLOTGEOMETRY = true
    mu = 3e10
    nu = 0.25
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

    plotbcsUeff(elsbox, bcsbox, Ueffbox, "Box")
    xgrid, ygrid = obsgrid(boxL+offset, boxB+offset, boxR-offset, boxT-offset, npts)
    Ubox, Sbox = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, 1:elsbox.endidx,
                                 Ueffbox[1:2:end], Ueffbox[2:2:end], mu, nu)
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ubox, Sbox, "Box")
    return

    
    # Element geometries and data structures for the layered box case
    elslayer = Elements(Int(1e5))
    nside = 40
    x1, y1, x2, y2 = discretizedline(-30e3, -15e3, 30e3, -15e3, nside) # Bottom
    addelsez!(elslayer, x1, y1, x2, y2, "B1")
    x1, y1, x2, y2 = discretizedline(30e3, -15e3, 30e3, 0, nside) # Right hand side
    addelsez!(elslayer, x1, y1, x2, y2, "R1")
    x1, y1, x2, y2 = discretizedline(30e3, 0, -30e3, 0, nside) # Top
    addelsez!(elslayer, x1, y1, x2, y2, "T1")
    x1, y1, x2, y2 = discretizedline(-30e3, 0, -30e3, -15e3, nside) # Left hand side
    addelsez!(elslayer, x1, y1, x2, y2, "L1")    
    x1, y1, x2, y2 = discretizedline(-30e3, -30e3, 30e3, -30e3, nside) # Bottom
    addelsez!(elslayer, x1, y1, x2, y2, "B2")
    x1, y1, x2, y2 = discretizedline(30e3, -30e3, 30e3, -15e3, nside) # Right hand side
    addelsez!(elslayer, x1, y1, x2, y2, "R2")
    x1, y1, x2, y2 = discretizedline(30e3, -15e3, -30e3, -15e3, nside) # Top
    addelsez!(elslayer, x1, y1, x2, y2, "T2")
    x1, y1, x2, y2 = discretizedline(-30e3, -15e3, -30e3, -30e3, nside) # Left hand side
    addelsez!(elslayer, x1, y1, x2, y2, "L2")
    idxlayer = getidxdict(elslayer)
    PLOTGEOMETRY &&  plotgeometry(elslayer, "Layer boundaries and normals")
    
    
    # Two layer box solutions
    TH = zeros(2*elslayer.endidx, 2*elslayer.endidx)
    C1idx = [idxlayer["F1"] ; idxlayer["B1"] ; idxlayer["R1"] ; idxlayer["T1"] ; idxlayer["L1"]]
    C2idx = [idxlayer["B2"] ; idxlayer["R2"] ; idxlayer["T2"] ; idxlayer["L2"]]
    
    # Boundaries shared by C1 and C2
    T_B1_C1, H_B1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["B1"], C1idx, mu, nu)
    T_T2_C2, H_T2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["T2"], C2idx, mu, nu)

    # C1 only boundaries
    T_F1_C1, T_H1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["F1"], C1idx, mu, nu)
    T_R1_C1, H_R1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["R1"], C1idx, mu, nu)
    T_T1_C1, H_T1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["T1"], C1idx, mu, nu)
    T_L1_C1, H_L1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["L1"], C1idx, mu, nu)

    # C2 only boundaries
    T_B2_C2, H_B2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["B2"], C2idx, mu, nu)
    T_R2_C2, H_R2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["R2"], C2idx, mu, nu)
    T_L2_C2, H_L2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["L2"], C2idx, mu, nu)

    # Place submatrices into larger matrix
    TH[1:80, 1:322] = T_B1_C1 # Displacements at B1 due to C1 (including F1)
    TH[1:80, 323:642] = -T_T2_C2 # Displacements at T2 due to C2
    TH[81:160, 1:322] = H_B1_C1 # Tractions at B1 due to C1 (including F1)
    TH[81:160, 323:642] = H_T2_C2 # Tractions at T2 due to C2
    TH[161:162, 1:322] = T_F1_C1 # Displacments at F1 due in C1 (including F1)
    TH[163:242, 1:322] = H_R1_C1 # Tractions at R1 due in C1 (including F1)
    TH[243:322, 1:322] = H_T1_C1 # Tractions at T1 due in C1 (including F1)
    TH[323:402, 1:322] = H_L1_C1 # Tractions at L1 due in C1 (including F1)
    TH[403:482, 323:642] = T_B2_C2 # Displacments at B2 due in C2
    TH[483:562, 323:642] = H_R2_C2 # Tractions at R2 due in C2
    TH[563:642, 323:642] = H_L2_C2 # Tractions at L2 due in C2
    bcslayer = zeros(2*elslayer.endidx)
    # bcslayer[161:162] .= 0.5
    bcslayer[1:2] .= 0.5

    @show cond(TH)

    matshow(log10.(abs.(TH)))
    
    figure()
    quiver(elslayer.xcenter[1:elslayer.endidx],
           elslayer.ycenter[1:elslayer.endidx],
           bcslayer[1:2:end], bcslayer[2:2:end])
    title("boundary conditions, WtF")
    
    # Solve BEM for effective displacements
    Uefflayer = inv(TH) \ bcslayer

    figure()
    quiver(elslayer.xcenter[1:elslayer.endidx],
           elslayer.ycenter[1:elslayer.endidx],
           Uefflayer[1:2:end], Uefflayer[2:2:end])
    
    return
    
    # xgrid, ygrid = obsgrid(-30e3+offset, -30e3+offset, 30e3-offset, -1-offset, npts)
    # xgrid1, ygrid1 = obsgrid(-30e3+offset, -15e3+offset, 30e3-offset, -1-offset, npts)
    # xgrid2, ygrid2 = obsgrid(-30e3+offset, -30e3+offset, 30e3-offset, -15e3-offset, npts)

    # Ulayer1, Slayer1 = constdispstress(slip2dispstress, xgrid1, ygrid1,
    #                                    elslayer, C1idx,
    #                                    Uefflayer[1:2:end], Uefflayer[2:2:end],
    #                                    mu, nu)
    # plotfields(elslayer, reshape(xgrid1, npts, npts), reshape(ygrid1, npts, npts),
    #            Ulayer1, Slayer1, "Layer 1 - no fault")

    # # Layer 2
    # Ulayer2, Slayer2 = constdispstress(slip2dispstress, xgrid2, ygrid2,
    #                                    elslayer, C2idx,
    #                                    Uefflayer[1:2:end], Uefflayer[2:2:end],
    #                                    mu, nu)
    # plotfields(elslayer, reshape(xgrid2, npts, npts), reshape(ygrid2, npts, npts),
    #            Ulayer2, Slayer2, "Layer 2 - no fault")

end
dislocationlayers()
