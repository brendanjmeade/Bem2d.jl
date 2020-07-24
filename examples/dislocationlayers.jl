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
    plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               U, S, "BEM halfspace")

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
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ubox, Sbox, "Box")


    #
    # Two layer box solutions
    #
    bcslayer = zeros(2*elslayer.endidx)
    bcslayer[4*nside:4*nside+1] .= 0.5
    figure()
    plot(bcslayer)
    show()
    
end
dislocationinabox()
