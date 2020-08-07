using Revise
using PyPlot
using PyCall
using Statistics
using LinearAlgebra
using Infiltrator
using Bem2d


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
    PLOTGEOMETRY = false
    mu = 3e10
    nu = 0.25
    
    # Element geometries and data structures for the box case
    elsbox = Elements(Int(1e5))
    nfault = 1
    boxbottom = -30000
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elsbox, x1, y1, x2, y2, "F")
    nside = 40
    x1, y1, x2, y2 = discretizedline(-30000, boxbottom, 30000, boxbottom, nside) # Bottom
    addelsez!(elsbox, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(30000, boxbottom, 30000, 0, nside) # Right hand side
    addelsez!(elsbox, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(30000, 0, -30000, 0, nside) # Top
    addelsez!(elsbox, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-30000, 0, -30000, boxbottom, nside) # Left hand side
    addelsez!(elsbox, x1, y1, x2, y2, "L")
    idxbox = getidxdict(elsbox)

    if PLOTGEOMETRY
        figure()
        for n in ["fault", "B", "R", "T", "L"]
            plot([elsbox.x1[idxbox[n]], elsbox.x2[idxbox[n]]  ],
                 [elsbox.y1[idxbox[n]], elsbox.y2[idxbox[n]]],
                 "-r")
            quiver(elsbox.xcenter[idxbox[n]],
                   elsbox.ycenter[idxbox[n]],
                   elsbox.xnormal[idxbox[n]],
                   elsbox.ynormal[idxbox[n]],
                   width=1e-3, scale = 1e2)
        end
        title("Box boundaries and normals")
        gca().set_aspect("equal")
    end
    
    # Element geometries and data structures for the layered box case
    elslayer = Elements(Int(1e5))
    nfault = 1
    nside = 40
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elslayer, x1, y1, x2, y2, "F1")
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

    if PLOTGEOMETRY
        scale = 500
        figure(figsize=(15,10))
        for i in 1:elslayer.endidx
            plot([elslayer.x1[i], elslayer.x2[i]],
                 [elslayer.y1[i], elslayer.y2[i]],
                 "-g")
            quiver(elslayer.xcenter[i],
                   elslayer.ycenter[i],
                   elslayer.xnormal[i],
                   elslayer.ynormal[i],
                   width=1e-3, scale=1e2)
            text(elslayer.xcenter[i] + scale * elslayer.xnormal[i],
                 elslayer.ycenter[i] + scale * elslayer.ynormal[i],
                 string(i), fontsize=5)
        end
        title("Layered")
        gca().set_aspect("equal")
    end
        
    # Box BEM problem
    T_FB_FBRTL, H_FB_FBRTL = PUTC(slip2dispstress, elsbox,
                                  [idxbox["F"]; idxbox["B"]],
                                  1:elsbox.endidx, mu, nu)
    T_RTL_FBRTL, H_RTL_FBRTL = PUTC(slip2dispstress, elsbox,
                                    [idxbox["R"]; idxbox["T"]; idxbox["L"]],
                                    1:elsbox.endidx, mu, nu)
    T_B_FBRTL, H_B_FBRTL = PUTC(slip2dispstress, elsbox, idxbox["B"],
                                      1:elsbox.endidx, mu, nu)
    bcsbox = zeros(2*elsbox.endidx)
    bcsbox[1:2] .= 0.5
    TH = [T_FB_FBRTL ; H_RTL_FBRTL]
    Ueffbox = TH \ bcsbox
    figure()
    quiver(elsbox.xcenter[1:elsbox.endidx],
           elsbox.ycenter[1:elsbox.endidx],
           Ueffbox[1:2:end], Ueffbox[2:2:end])
    title("Ueff for the total box")
    
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
    TH[81:160, 323:642] = H_T2_C2 # Tractions at T2 due to C2
    TH[81:160, 1:322] = H_B1_C1 # Tractions at B1 due to C1 (including F1)
    TH[161:162, 1:322] = T_F1_C1 # Displacments at F1 due in C1 (including F1)
    TH[163:242, 1:322] = H_R1_C1 # Tractions at R1 due in C1 (including F1)
    TH[243:322, 1:322] = H_T1_C1 # Tractions at T1 due in C1 (including F1)
    TH[323:402, 1:322] = H_L1_C1 # Tractions at L1 due in C1 (including F1)
    TH[403:482, 323:642] = T_B2_C2 # Displacments at B2 due in C2
    TH[483:562, 323:642] = H_R2_C2 # Tractions at R2 due in C2
    TH[563:642, 323:642] = H_L2_C2 # Tractions at L2 due in C2
    bcslayer = zeros(2*elslayer.endidx)
    bcslayer[161] = 0.5
    bcslayer[162] = 0.5
    @show cond(TH)

    # # Let's do the upper layer only. Should be able to compare with the box solution
    # TH = zeros(2 * length(C1idx), 2 * length(C1idx))

    # @show size(TH)
    # TH = [T_FB_FBRTL ; H_RTL_FBRTL] # This is a good reference set of partials from the box problem
    # TH[1:2, :] = T_F1_C1
    # TH[3:82, :] = T_B1_C1 # Displacments at B1 due in C1 (including F1)
    # TH[83:162, :] = H_R1_C1 # Tractions at R1 due in C1 (including F1)
    # TH[163:242, :] = H_T1_C1 # Tractions at T1 due in C1 (including F1)
    # TH[243:322, :] = H_L1_C1 # Tractions at L1 due in C1 (including F1)
    # bcsul = zeros(2*length(C1idx))
    # bcsul[1] = 0.5
    # bcsul[2] = 0.5
    # Uefful = TH \ bcsul

    # figure()
    # quiver(elslayer.xcenter[C1idx], elslayer.ycenter[C1idx],
    #        Uefful[1:2:end], Uefful[2:2:end])
    # title("Ueff for upper layer only")

    # return

    
    bcstemp = zeros(2*length(C1idx))
    bcstemp[161] = 0.5
    bcstemp[162] = 0.5
    
    figure()
    quiver(elslayer.xcenter[1:elslayer.endidx],
           elslayer.ycenter[1:elslayer.endidx],
           bcslayer[1:2:end], bcslayer[2:2:end])
    
    # Solve BEM for effective displacements
    Uefflayer = inv(TH) \ bcslayer
    matshow(log10.(abs.(TH))); colorbar()

    figure()
    quiver(elslayer.xcenter[1:elslayer.endidx],
           elslayer.ycenter[1:elslayer.endidx],
           Uefflayer[1:2:end], Uefflayer[2:2:end])


    # Try a little demo proble in the uper layer only
    
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
