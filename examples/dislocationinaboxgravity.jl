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
    lambda = 30e9
    nu = 0.25
    g = 9.81
    rho = 2700
    alpha = 1e-7 # scalar preconditioner
    npts = 50
    offset = 1
    xgrid, ygrid = obsgrid(-30e3+offset, -20e3+offset, 30e3-offset, -1-offset, npts)

    # Element geometries and data structures for the box case
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

    # Box BEM problem (dislocation only, no gravity)
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
    THmat = [T_faultB_faultBRTL ; alpha .* H_RTL_faultBRTL]
    Ueffbox = inv(THmat) * bcsbox
    Ubox, Sbox = constdispstress(slip2dispstress, xgrid, ygrid, elsbox,
        collect(1:1:elsbox.endidx), Ueffbox[1:2:end], Ueffbox[2:2:end], mu, nu)
    # plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
    #            Ubox, Sbox, "Box")


    #
    # Box BEM problem (no dislocation, gravity only)
    #
    Ug, Sg = gravityparticularfunctions(elsbox.xcenter[idxbox["B"]],
                                        elsbox.ycenter[idxbox["B"]], g, rho, lambda, mu)    
    bcsboxgravity = zeros(2 * elsbox.endidx - 2)
    bcsboxgravity[1:2:2*nside] = Ug[:, 1]
    bcsboxgravity[2:2:2*nside] = Ug[:, 2]
    bcsboxgravity *= -1
    Ueffboxparticular = inv(THmat[3:end, 3:end]) * bcsboxgravity # Clip off the fault part    
    Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, collect(2:1:elsbox.endidx),
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)
    Uint, Sint = gravityparticularfunctions(xgrid, ygrid, g, rho, lambda, mu)
    U = @. Ucomp + Uint
    S = @. Scomp + Sint
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts), U, S, "Complementary + Particular solutions")

    
    #
    # Box BEM problem (dislocation and gravity)
    #

    # Particular solution and effective boundary conditions
    Ug, Sg = gravityparticularfunctions(elsbox.xcenter[1:1:elsbox.endidx],
                                        elsbox.ycenter[1:1:elsbox.endidx], g, rho, lambda, mu)    
    # bcsboxgravity = zeros(2 * length([idxbox["fault"]; idxbox["B"]]))
    bcsboxgravity = zeros(2 * elsbox.endidx)
    @infiltrate
    
    # bcsboxgravity[1:2:6*nels] = UB[:, 1] # Bottom boundary (x-component)
    # bcsboxgravity[2:2:6*nels] = UB[:, 2] # Bottom boundary (y-component)
    # bcsboxgravity *= -1 # This is neccesary for the right answer and is consistent with derivation

    # BEM solve to get particular solution
    # Ueffparticular = inv(THmat) * bcsgravity

    # # Evaluate and plot interior solution
    # Uinteriorcomplementary, Sinteriorcomplementary = quaddispstress(slip2dispstress, x, y, els, bcidxall, quadstack(Ueffparticular[1:2:end]), quadstack(Ueffparticular[2:2:end]), mu, nu)
    # Uinteriorparticular, Sinteriorparticular = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    # U = @. Uinteriorcomplementary + Uinteriorparticular
    # S = @. Sinteriorcomplementary + Sinteriorparticular
    # # plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts), U, S, "Complementary + Particular solutions")


end
dislocationinabox()
