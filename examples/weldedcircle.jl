using Revise
using PyPlot
using PyCall
using Statistics
using LinearAlgebra
using Infiltrator
using IterativeSolvers
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
    weldedcircle()

Comparing homogeneous and welded circle BEM solutions
"""
function weldedcircle()
    close("all")
    PLOTGEOMETRY = true
    # mu = 3e10
    # nu = 0.25
    mu1 = 3e10
    nu1 = 0.25
    mu2 = 1.0 * mu1
    nu2 = nu1

    npts = 30
    offset = 100 # meters
    
    # Element geometries and data structures for the homogeneous circle case
    els1 = Elements(Int(1e5))
    nels = 20
    r = 10e3
    x1, y1, x2, y2 = discretizedarc(deg2rad(0), deg2rad(180), r, nels)
    addelsez!(els1, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedarc(deg2rad(180), deg2rad(360), r, nels)
    addelsez!(els1, x1, y1, x2, y2, "B")
    idx1 = getidxdict(els1)
    # PLOTGEOMETRY && plotgeometry(els1, "Circle boundaries and normals")

    T_B_BT, _ = PUTC(slip2dispstress, els1, idx1["B"], 1:els1.endidx, mu1, nu1)
    _, H_T_BT = PUTC(slip2dispstress, els1, idx1["T"], 1:els1.endidx, mu1, nu1)
    bcs1 = zeros(2*els1.endidx)
    bcs1[2*10] = 1.0
    bcs1[2*11] = 1.0    
    Ueff1 = [H_T_BT ; T_B_BT] \ bcs1
    # plotbcsUeff(els1, bcs1, Ueff1, "homogeneous circle")
    xgrid1, ygrid1 = obsgrid(-r, -r, r, r, npts)
    U1, S1 = constdispstress(slip2dispstress, xgrid1, ygrid1, els1, 1:els1.endidx,
                             Ueff1[1:2:end], Ueff1[2:2:end], mu1, nu1)
    # plotfields(els1, reshape(xgrid1, npts, npts), reshape(ygrid1, npts, npts),
    #            U1, S1, "homogeneous circle")
  
    # Element geometries and data structures for the welded circle case
    els2 = Elements(Int(1e5))
    nels = 20
    r = 10e3
    x1, y1, x2, y2 = discretizedarc(deg2rad(0), deg2rad(180), r, nels)
    addelsez!(els2, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-r, 0, r, 0, nels)
    addelsez!(els2, x1, y1, x2, y2, "midT")
    x1, y1, x2, y2 = discretizedarc(deg2rad(180), deg2rad(360), r, nels)
    addelsez!(els2, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(-r, 0, r, 0, nels) # Peculiar
    addelsez!(els2, x2, y2, x1, y1, "midB")
    idx2 = getidxdict(els2)
    # PLOTGEOMETRY && plotgeometry(els2, "Welded circle boundaries and normals")
    
    TH = zeros(2 * els2.endidx, 2 * els2.endidx)
    T_midT_TmidT, H_midT_TmidT = PUTC(slip2dispstress, els2, idx2["midT"],
                                      [idx2["T"]; idx2["midT"]],
                                      mu1, nu1)
    T_midB_BmidB, H_midB_BmidB = PUTC(slip2dispstress, els2,
                                      idx2["midB"], [idx2["B"]; idx2["midB"]],
                                      mu2, nu2)
    _, H_T_TmidT = PUTC(slip2dispstress, els2, idx2["T"],
                        [idx2["T"]; idx2["midT"]], mu1, nu1)
    T_B_BmidB, _ = PUTC(slip2dispstress, els2, idx2["B"],
                        [idx2["B"]; idx2["midB"]], mu2, nu2)

    TH[1:40, 1:80] = T_midT_TmidT
    TH[1:40, 81:160] = -T_midB_BmidB
    TH[41:80, 1:80] = H_midT_TmidT
    TH[41:80, 81:160] = H_midB_BmidB
    TH[81:120, 1:80] = H_T_TmidT
    TH[121:160, 81:160] = T_B_BmidB

    matshow(log10.(abs.(TH)))
    title("condition number = " * string(cond(TH)))
    colorbar()
    
    bcs2 = zeros(2*els2.endidx)
    # bcs2[2*10] = 1.0 # These are the initial indices
    # bcs2[2*11] = 1.0
    bcs2[2*50] = 1.0 # These are the indices after welding
    bcs2[2*51] = 1.0

    # Direct solve
    Ueff2 = TH \ bcs2
    
    # Try iterative solver
    # Ueff2gmres, history = gmres(TH, bcs2, log=true, verbose=true)
    # Ueff2idrs, history = idrs(TH, bcs2; s = 8, log=true, verbose=true)
    # Ueff2 = Ueff2idrs
    # figure()
    # plot(Ueff2, "xb")
    # plot(Ueff2gmres, "+r")
    UeffT2 = Ueff2[1:80]
    UeffB2 = Ueff2[81:160]
    # plotbcsUeff(els2, bcs2, Ueff2, "welded circle")
    
    
    # Internal evaluation
    Tidx = findall(x -> x > 0, ygrid1)
    Bidx = findall(x -> x < 0, ygrid1)
    Tx = xgrid1[Tidx]
    Ty = ygrid1[Tidx]
    Bx = xgrid1[Bidx]
    By = ygrid1[Bidx]
    UT2, ST2 = constdispstress(slip2dispstress, Tx, Ty, els2,
                               [idx2["T"]; idx2["midT"]],
                               UeffT2[1:2:end], UeffT2[2:2:end], mu1, nu1)
    UB2, SB2 = constdispstress(slip2dispstress, Bx, By, els2,
                               [idx2["B"]; idx2["midB"]],
                               UeffB2[1:2:end], UeffB2[2:2:end], mu2, nu2)
    
    # xgridT2, ygridT2 = obsgrid(-r, 344.827, r, r, npts)
    # UT2, ST2 = constdispstress(slip2dispstress, xgridT2, ygridT2, els2,
    #                            [idx2["T"]; idx2["midT"]],
    #                            UeffT2[1:2:end], UeffT2[2:2:end], mu, nu)
    # # plotfields(els2, reshape(xgridT2, npts, npts), reshape(ygridT2, npts, npts),
    # #            UT2, ST2, "homogeneous circle (top)")

    # xgridB2, ygridB2 = obsgrid(-r, -r, r, -344.827, npts)
    # UB2, SB2 = constdispstress(slip2dispstress, xgridB2, ygridB2, els2,
    #                            [idx2["B"]; idx2["midB"]],
    #                            UeffB2[1:2:end], UeffB2[2:2:end], mu, nu)
    # plotfields(els2, reshape(xgridB2, npts, npts), reshape(ygridB2, npts, npts),
    #            UB2, SB2, "homogeneous circle (bottom)")

    # Quick and dirty reshaping
    figure(figsize=(20, 10))
    subplot(1, 2, 1)
    quiver(xgrid1[:], ygrid1[:], U1[:, 1], U1[:, 2],
           width=1e-3, scale=1e-6)
    gca().set_aspect("equal")

    subplot(1, 2, 2)
    quiver(Tx[:], Ty[:], UT2[:, 1], UT2[:, 2],
           width=1e-3, scale=1e-6)
    quiver(Bx[:], By[:], UB2[:, 1], UB2[:, 2],
           width=1e-3, scale=1e-6)    
    gca().set_aspect("equal")

end
weldedcircle()
