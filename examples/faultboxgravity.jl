using Revise
using PyPlot
using PyCall
using Statistics
using LinearAlgebra
using Infiltrator
using Bem2d


"""
    faultboxgravity()

Comparing half-space and dislocaiton in a box solutions
"""
function faultboxgravity()
    close("all")
    mu = 30e9
    lambda = 30e9
    nu = 0.25
    g = 9.81
    rho = 2700
    alpha = 1e-7 # scalar preconditioner
    npts = 100
    offset = 1
    B = -15e3 # Bottom
    R = 15e3 # Right
    T = 0e3 # Top
    L = -15e3 # Left
    xgrid, ygrid = obsgrid(L+offset, B+offset, R-offset, T-offset, npts)

    # Element geometries and data structures for the box case
    elsbox = Elements(Int(1e5))
    nfault = 1
    nside = 20
    x1, y1, x2, y2 = discretizedline(L, B, R, B, nside) # Bottom
    addelsez!(elsbox, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(R, B, R, T, nside) # Right hand side
    addelsez!(elsbox, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(R, T, L, T, nside) # Top
    addelsez!(elsbox, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(L, T, L, B, nside) # Left hand side
    addelsez!(elsbox, x1, y1, x2, y2, "L")
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elsbox, x1, y1, x2, y2, "F")
    idxbox = getidxdict(elsbox)

    # Fault contribution
    T_TB_F, H_TB_F = PUTC(slip2dispstress, elsbox,
                          [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                          idxbox["F"], mu, nu)
    Fslip = [10; 10]; # y-direction slip only
    Uslip = T_TB_F * Fslip
    Tslip = H_TB_F * Fslip

    # Kernels and solve
    T_B_BRTL, H_B_BRTL = PUTC(slip2dispstress, elsbox, idxbox["B"],
                       [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                       mu, nu)
    T_R_BRTL, H_R_BRTL = PUTC(slip2dispstress, elsbox, idxbox["R"],
                       [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                       mu, nu)
    T_L_BRTL, H_L_BRTL = PUTC(slip2dispstress, elsbox, idxbox["L"],
                       [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                       mu, nu)
    _, H_T_BRTL = PUTC(slip2dispstress, elsbox, idxbox["T"],
                       [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                       mu, nu)

    bcsbox = zeros(8 * nside)
    bcsbox[1:2*nside] = -Uslip[1:2*nside] # Bottom
    bcsbox[2*nside+1:4*nside] = -Tslip[2*nside+1:4*nside] # Right
    bcsbox[4*nside+1:6*nside] = -Tslip[4*nside+1:6*nside] # Top
    bcsbox[6*nside+1:8*nside] = -Tslip[6*nside+1:8*nside] # Left
    THbox = [T_B_BRTL ; H_R_BRTL; H_T_BRTL ; H_L_BRTL]
    Ueffbox = THbox \ bcsbox

    # Volume evaluation
    UTB, STB = constdispstress(slip2dispstress, xgrid, ygrid, elsbox,
                               [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                               Ueffbox[1:2:end], Ueffbox[2:2:end], mu, nu)
    UF, SF = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, idxbox["F"],
                             Fslip[1:2:end], Fslip[2:2:end], mu, nu)
    Ufaultonly = UTB .+ UF
    Sfaultonly = STB .+ SF
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ufaultonly, Sfaultonly, "(Fault only)")

    

    
    ###
    ### Box BEM problem (no dislocation, gravity only)
    ###
    Ug, Sg = gravityparticularfunctions(elsbox.xcenter[idxbox["B"]],
                                        elsbox.ycenter[idxbox["B"]],
                                        g, rho, lambda, mu)    
    bcsboxgravity = zeros(2 * elsbox.endidx - 2)
    bcsboxgravity[1:2:2*nside] = Ug[:, 1]
    bcsboxgravity[2:2:2*nside] = Ug[:, 2]
    bcsboxgravity *= -1
    Ueffboxparticular = THbox \ bcsboxgravity # Clip off the fault part    

    Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox,
                                   [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)
    Uint, Sint = gravityparticularfunctions(xgrid, ygrid, g, rho, lambda, mu)
    Ugravityonly = @. Ucomp + Uint
    Sgravityonly = @. Scomp + Sint
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ugravityonly, Sgravityonly, "(Gravity only)")



    
    ###
    ### Box BEM problem (dislocation and gravity)
    ###
    bcscombined = bcsboxgravity .+ bcsbox    
    Ueffcombined = THbox \ bcsboxgravity # Clip off the fault part    
    Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox,
                                   [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)

    Ugravityandfault = Ucomp .+ Uint .+ UF
    Sgravityandfault = Scomp .+ Sint .+ SF
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ugravityandfault, Sgravityandfault, "(Gravity + fault)")
    # Difference between gravity alone and gravity + fault
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ugravityonly.-Ugravityandfault, Sgravityonly.-Sgravityandfault, "(Gravity only) - (Gravity + fault)")
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               UF, SF, "Fullspace fault")
end
faultboxgravity()
