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
    npts = 50
    offset = 1
    xgrid, ygrid = obsgrid(-30e3+offset, -20e3+offset, 30e3-offset, -1-offset, npts)

    # Element geometries and data structures for the box case
    elsbox = Elements(Int(1e5))
    nfault = 1
    nside = 20
    x1, y1, x2, y2 = discretizedline(-30000, -20000, 30000, -20000, nside) # Bottom
    addelsez!(elsbox, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(30000, -20000, 30000, 0, nside) # Right hand side
    addelsez!(elsbox, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(30000, 0, -30000, 0, nside) # Top
    addelsez!(elsbox, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-30000, 0, -30000, -20000, nside) # Left hand side
    addelsez!(elsbox, x1, y1, x2, y2, "L")
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elsbox, x1, y1, x2, y2, "F")
    idxbox = getidxdict(elsbox)

    # Fault contribution
    T_TB_F, H_TB_F = PUTC(slip2dispstress, elsbox,
                          [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                          idxbox["F"], mu, nu)
    Fslip = [1; 1]; # y-direction slip only
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
    U1 = UTB .+ UF
    S1 = STB .+ SF
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               U1, S1, "(Fault only)")
    
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
    U = @. Ucomp + Uint
    S = @. Scomp + Sint
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               U, S, "(Gravity only)")

    ###
    ### Box BEM problem (dislocation and gravity)
    ###
    bcscombined = bcsboxgravity .+ bcsbox    
    Ueffcombined = THbox \ bcsboxgravity # Clip off the fault part    
    Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox,
                                   [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)

    Ucombined = Ucomp .+ Uint .+ UF
    Scombined = Scomp .+ Sint .+ SF
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ucombined, Scombined, "(Gravity + fault)")

    # Difference between gravity alone and gravity + fault
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ucombined.-U, Scombined.-S, "(Gravity only) - (Gravity + fault)")

    # Difference between gravity alone and gravity + fault and fault alone
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ucombined.-U.-U1, Scombined.-S.-S1, "(Gravity only) - (Gravity + fault) - (fault only)")

    
    ###
    ### Box BEM problem (dislocation and gravity)
    ###
    # Ug, Sg = gravityparticularfunctions(elsbox.xcenter[1:1:elsbox.endidx],
    #                                     elsbox.ycenter[1:1:elsbox.endidx], g, rho, lambda, mu)    
    # bcsboxgravity = zeros(2 * elsbox.endidx)
    # bcsboxgravity[1:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 1]
    # bcsboxgravity[2:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 2]
    # bcsboxgravity *= -1
    # bcsboxgravity[1:2] .= 0.0 # Eliminate gravity here
    # Ueffboxparticular = inv(THmat) * bcsboxgravity    
    # Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, collect(1:1:elsbox.endidx),
    #                                Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)
    # Uint, Sint = gravityparticularfunctions(xgrid, ygrid, g, rho, lambda, mu)
    # Udg = @. Ucomp + Uint
    # Sdg = @. Scomp + Sint
    # plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
    #            Udg, Sdg, "Gravity on box edges and dislocation")
    
    # #
    # # Box BEM gravity with a dislocation
    # #
    # bcsboxgravity = zeros(2 * elsbox.endidx)
    # bcsboxgravity[1:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 1]
    # bcsboxgravity[2:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 2]
    # bcsboxgravity *= -1
    # # bcsboxgravity[1:2] .+= 0.5
    # bcsboxgravity[1:2] .= 0.5 # Add fault slip and eliminate gravity here
    # @show cond(THmat)
    # Ueffboxparticular = inv(THmat) * bcsboxgravity    
    # Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, collect(1:1:elsbox.endidx),
    #                                Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)
    # Uint, Sint = gravityparticularfunctions(xgrid, ygrid, g, rho, lambda, mu)
    # Udg2 = @. Ucomp + Uint
    # Sdg2 = @. Scomp + Sint
    # plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
    #            Udg2, Sdg2, "Gravity on box edges and dislocation")
end
faultboxgravity()
