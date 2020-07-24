using Revise
using PyPlot
using PyCall
using Statistics
using LinearAlgebra
using Infiltrator
using Bem2d


"""
    dislocationinabox()

Comparing half-space and dislocaiton in a box solutions
"""
function dislocationinabox()
    close("all")
    mu = 30e9
    nu = 0.25
    alpha = 1e-7 # scalar preconditioner
    npts = 50
    offset = 1
    xgrid, ygrid = obsgrid(-30e3+offset, -20e3+offset, 30e3-offset, -1-offset, npts)

    #
    # Element geometries and data structures for half space approximation
    #
    elshalfspace = Elements(Int(1e5))
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elshalfspace, x1, y1, x2, y2, "fault")
    nsurf = 100
    x1, y1, x2, y2 = discretizedline(-100e3, 0, 100e3, 0, nsurf) # Free surface
    addelsez!(elshalfspace, x1, y1, x2, y2, "surf")
    idxhalfspace = getidxdict(elshalfspace)

    #
    # Element geometries and data structures for the box case
    #
    elsbox = Elements(Int(1e5))
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elsbox, x1, y1, x2, y2, "fault")
    nside = 40
    x1, y1, x2, y2 = discretizedline(-30000, -20000, 30000, -20000, nside) # Bottom
    addelsez!(elsbox, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(30000, -20000, 30000, 0, nside) # Right hand side
    addelsez!(elsbox, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(30000, 0, -30000, 0, nside) # Top
    addelsez!(elsbox, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-30000, 0, -30000, -20000, nside) # Left hand side
    addelsez!(elsbox, x1, y1, x2, y2, "L")
    idxbox = getidxdict(elsbox)

    #
    # Element geometries and data structures for the layered box case
    #
    elslayer = Elements(Int(1e5))
    nfault = 1
    nside = 40
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elslayer, x1, y1, x2, y2, "F1")
    x1, y1, x2, y2 = discretizedline(-30000, -15000, 30000, -15000, nside) # Bottom
    addelsez!(elslayer, x1, y1, x2, y2, "B1")
    x1, y1, x2, y2 = discretizedline(30000, -15e3, 30000, 0, nside) # Right hand side
    addelsez!(elslayer, x1, y1, x2, y2, "R1")
    x1, y1, x2, y2 = discretizedline(30000, 0, -30000, 0, nside) # Top
    addelsez!(elslayer, x1, y1, x2, y2, "T1")
    x1, y1, x2, y2 = discretizedline(-30000, 0, -30000, -15000, nside) # Left hand side
    addelsez!(elslayer, x1, y1, x2, y2, "L1")
    x1, y1, x2, y2 = discretizedline(30000, -20e3, -30000, -20e3, nside) # Top
    addelsez!(elslayer, x1, y1, x2, y2, "B2")
    x1, y1, x2, y2 = discretizedline(30000, -15e3, 30000, -20e3, nside) # Right hand side
    addelsez!(elslayer, x1, y1, x2, y2, "R2")
    x1, y1, x2, y2 = discretizedline(-30000, -15000, 30000, -15000, nside) # Bottom
    addelsez!(elslayer, x1, y1, x2, y2, "T2")
    x1, y1, x2, y2 = discretizedline(-30000, -15e3, -30000, -20e3, nside) # Left hand side
    addelsez!(elslayer, x1, y1, x2, y2, "L2")
    idxlayer = getidxdict(elslayer)
    plotelements(elslayer)

    
    #
    # Halfspace BEM solution
    #
    T_fault_all, _ = PUTC(slip2dispstress, elshalfspace, idxhalfspace["fault"],
                          collect(1:1:elshalfspace.endidx), mu, nu)
    _, H_surf_all = PUTC(slip2dispstress, elshalfspace, idxhalfspace["surf"],
                         collect(1:1:elshalfspace.endidx), mu, nu)
    bcs = zeros(2 * elshalfspace.endidx)
    bcs[1:2*nfault] .= 0.5
    Ueff = inv([T_fault_all; alpha .* H_surf_all]) * bcs
    U, S = constdispstress(slip2dispstress, xgrid, ygrid, elshalfspace,
                           collect(1:1:elshalfspace.endidx),
                           Ueff[1:2:end], Ueff[2:2:end], mu, nu)
    # plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
    #            U, S, "BEM halfspace")

    #
    # Box BEM problem
    #
    T_faultB_faultBRTL, H_pfaultB_faultBRTL = PUTC(slip2dispstress, elsbox,
                                                   [idxbox["fault"]; idxbox["B"]],
                                                   collect(1:1:elsbox.endidx),
                                                   mu, nu)
    T_RTL_faultBRTL, H_RTL_faultBRTL = PUTC(slip2dispstress, elsbox,
                                            [idxbox["R"]; idxbox["T"]; idxbox["L"]],
                                            collect(1:1:elsbox.endidx),
                                            mu, nu)
    bcsbox = zeros(2*elsbox.endidx)
    bcsbox[1:2] .= 0.5
    Ueffbox = inv([T_faultB_faultBRTL ; alpha .* H_RTL_faultBRTL]) * bcsbox
    Ubox, Sbox = constdispstress(slip2dispstress, xgrid, ygrid, elsbox,
        collect(1:1:elsbox.endidx), Ueffbox[1:2:end], Ueffbox[2:2:end], mu, nu)
    # plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
    #            Ubox, Sbox, "Box")


    # #
    # # Two layer box solutions
    # #
    # TH = zeros(2*elslayer.endidx, 2*elslayer.endidx)
    # C1idx = [idxlayer["F1"] ; idxlayer["B1"] ; idxlayer["R1"] ; idxlayer["T1"] ; idxlayer["L1"]]
    # C2idx = [idxlayer["B2"] ; idxlayer["R2"] ; idxlayer["T2"] ; idxlayer["L2"]]
    # C1idxlong = collect(minimum(C1idx):1:2*maximum(C1idx))
    # C2idxlong = collect(2*minimum(C2idx):1:2*maximum(C2idx))
    
    # # Boundaries shared by C1 and C2
    # T_B1_C1, H_B1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["B1"], C1idx, mu, nu)
    # T_T2_C2, H_T2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["T2"], C2idx, mu, nu)

    # # C1 only boundaries
    # T_F1_C1, _ = PUTC(slip2dispstress, elslayer, idxlayer["F1"], C1idx, mu, nu)
    # _, H_R1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["R1"], C1idx, mu, nu)
    # _, H_T1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["T1"], C1idx, mu, nu)
    # _, H_L1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["L1"], C1idx, mu, nu)

    # #C2 only boundaries
    # T_B2_C2, _ = PUTC(slip2dispstress, elslayer, idxlayer["B2"], C2idx, mu, nu)
    # _, H_R2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["R2"], C2idx, mu, nu)
    # _, H_L2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["L2"], C2idx, mu, nu)

    # # Place submatrices into large matrix
    # TH[1:80, 1:322] = T_B1_C1
    # TH[81:160, 1:322] = alpha .* H_B1_C1
    # TH[1:80, 323:642] = T_T2_C2
    # TH[81:160, 323:642] = -alpha .* H_T2_C2
    # TH[161:162, 1:322] = T_F1_C1 # Fault
    # TH[163:242, 1:322] = alpha .* H_R1_C1
    # TH[243:322, 1:322] = alpha .* H_T1_C1
    # TH[323:402, 1:322] = alpha .* H_L1_C1
    # TH[403:482, 323:642] = T_B2_C2
    # TH[483:562, 323:642] = alpha .* H_R2_C2
    # TH[563:642, 323:642] = alpha .* H_L2_C2

    # bcslayer = zeros(2*elslayer.endidx)
    # bcslayer[1:2] .= 0.5

    # Uefflayer = inv(TH) * bcslayer
    # matshow(log10.(abs.(TH)))
    # colorbar()
    # title("TH")

    # matshow(log10.(abs.(inv(TH))))
    # colorbar()
    # title("TH")

    # figure()
    # quiver(elslayer.xcenter[1:elslayer.endidx], elslayer.ycenter[1:elslayer.endidx],
    #        Uefflayer[1:2:end], Uefflayer[2:2:end])


    # Trying without fault (condition number is all that matters)
    TH = zeros(2*elslayer.endidx-2, 2*elslayer.endidx-2)
    C1idx = [idxlayer["B1"] ; idxlayer["R1"] ; idxlayer["T1"] ; idxlayer["L1"]]
    C2idx = [idxlayer["B2"] ; idxlayer["R2"] ; idxlayer["T2"] ; idxlayer["L2"]]

    # Boundaries shared by C1 and C2
    T_B1_C1, H_B1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["B1"], C1idx, mu, nu)
    T_T2_C2, H_T2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["T2"], C2idx, mu, nu)

    # C1 only boundaries
    _, H_R1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["R1"], C1idx, mu, nu)
    _, H_T1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["T1"], C1idx, mu, nu)
    _, H_L1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["L1"], C1idx, mu, nu)

    #C2 only boundaries
    T_B2_C2, _ = PUTC(slip2dispstress, elslayer, idxlayer["B2"], C2idx, mu, nu)
    _, H_R2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["R2"], C2idx, mu, nu)
    _, H_L2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["L2"], C2idx, mu, nu)

    # Place submatrices into large matrix
    TH[1:80, 1:320] = T_B1_C1
    TH[81:160, 1:320] = -alpha .* H_B1_C1
    TH[1:80, 321:640] = T_T2_C2
    TH[81:160, 321:640] = alpha .* H_T2_C2
    TH[161:240, 1:320] = alpha .* H_R1_C1
    TH[241:320, 1:320] = alpha .* H_T1_C1
    TH[321:400, 1:320] = alpha .* H_L1_C1
    TH[401:480, 321:640] = T_B2_C2
    TH[481:560, 321:640] = alpha .* H_R2_C2
    TH[561:640, 321:640] = alpha .* H_L2_C2


    bcslayer = zeros(2*elslayer.endidx-2)
    bcslayer[281:282] .= 0.5

    Uefflayer = inv(TH) * bcslayer
    matshow(log10.(abs.(TH)))
    colorbar()
    title("TH")

    matshow(log10.(abs.(inv(TH))))
    colorbar()
    title("TH")

figure()
plot(log10.(abs.(Uefflayer)))
    figure()
quiver(elslayer.xcenter[2:elslayer.endidx],
       elslayer.ycenter[2:elslayer.endidx],
       Uefflayer[1:2:end], Uefflayer[2:2:end])

    # Try a displacement only bottom layer problem
    T_B2_C2, H_B2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["B2"], C2idx, mu, nu)
    H_R2_C2, H_R2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["R2"], C2idx, mu, nu)
    T_T2_C2, H_T2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["T2"], C2idx, mu, nu)
    H_L2_C2, H_L2_C2 = PUTC(slip2dispstress, elslayer, idxlayer["L2"], C2idx, mu, nu)
    TH = zeros(320, 320)
    TH[1:80, :] = T_B2_C2
    TH[81:160, :] = alpha .* H_R2_C2
    TH[161:240, :] = alpha .* H_T2_C2
    TH[241:320, :] = alpha .* H_L2_C2

    @show cond(TH)
    @show rank(TH)
    @show det(TH)

    # Try a displacement only top layer problem
    T_B1_C1, H_B1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["B1"], C1idx, mu, nu)
    H_R1_C1, H_R1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["R1"], C1idx, mu, nu)
    T_T1_C1, H_T1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["T1"], C1idx, mu, nu)
    H_L1_C1, H_L1_C1 = PUTC(slip2dispstress, elslayer, idxlayer["L1"], C1idx, mu, nu)
    TH = zeros(320, 320)
    TH[1:80, :] = T_B1_C1
    TH[81:160, :] = alpha .* H_R1_C1
    TH[161:240, :] = alpha .* H_T1_C1
    TH[241:320, :] = alpha .* H_L1_C1

    @show cond(TH)
    @show rank(TH)
    @show det(TH)


@infiltrate

end
dislocationinabox()