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
    B = -20e3 # Bottom
    R = 20e3 # Right
    T = 0e3 # Top
    L = -20e3 # Left

    # Element geometries and data structures for the box case
    elsbox = Elements(Int(1e5))
    nfault = 1
    nside = 20

    # From thrust fault example
    x1T, y1T, x2T, y2T = discretizedline(-20e3, 0, 20e3, 0, nside)
    y1T = -1e3 * atan.(x1T / 1e3)
    y2T = -1e3 * atan.(x2T / 1e3)
    xgrid, ygrid = obsgrid(L+offset, B+offset, R-offset, maximum(y1T)-offset, npts)

    x1, y1, x2, y2 = discretizedline(L, B, R, B, nside) # Bottom
    addelsez!(elsbox, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(R, B, R, minimum(y2T), nside) # Right hand side
    addelsez!(elsbox, x1, y1, x2, y2, "R")
    addelsez!(elsbox, x1T, y1T, x2T, y2T, "T")
    x1, y1, x2, y2 = discretizedline(L, maximum(y1T), L, B, nside) # Left hand side
    addelsez!(elsbox, x1, y1, x2, y2, "L")
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, -5e3, -5e3, nfault) # 45 degree dipping fault
    addelsez!(elsbox, x1, y1, x2, y2, "F")
    idxbox = getidxdict(elsbox)

    ###
    ### Fault only
    ### R, T, L are traction free bcs
    ### B is zero slip BC
    ### 
    T_TB_F, H_TB_F = PUTC(slip2dispstress, elsbox,
                          [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                          idxbox["F"], mu, nu)
    Fslip = [100; 100]; # y-direction slip only
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


    figure()
    field = sqrt.(Ufaultonly[:, 1].^2 + Ufaultonly[:, 2].^2)
    ncontours = 40
    xlim = [-20000 20000]
    ylim = [-20000 2000]
    scale = 1.0
    fieldmax = maximum(@.abs(field))
    contourf(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
             reshape(field, npts, npts), ncontours,
             vmin=-scale * fieldmax, vmax=scale * fieldmax,
             cmap=rycroftcmap())
    clim(-scale * fieldmax, scale * fieldmax)
    colorbar(fraction=0.020, pad=0.05, extend="both")
    contour(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
            reshape(field, npts, npts), ncontours,
            vmin=-scale * fieldmax, vmax=scale * fieldmax,
            linewidths=0.25, colors="k")
    title(title)
    plotelements(elsbox)

    gca().set_aspect("equal")
    gca().set_xlim([xlim[1], xlim[2]])
    gca().set_ylim([ylim[1], ylim[2]])
    # gca().set_xticks([xlim[1], xlim[2]])
    # gca().set_yticks([ylim[1], ylim[2]])
   
    return

    
    ###
    ### Box BEM problem (no dislocation, gravity only)
    ###
    # Displacement BCs for bottom
    Ug, Sg = gravityparticularfunctions(elsbox.xcenter[idxbox["B"]],
                                        elsbox.ycenter[idxbox["B"]],
                                        g, rho, lambda, mu)    
    bcsboxgravity = zeros(2 * elsbox.endidx - 2)
    bcsboxgravity[1:2:2*nside] = Ug[:, 1]
    bcsboxgravity[2:2:2*nside] = Ug[:, 2]
    
    # Traction BCs for top
    Ug, Sg = gravityparticularfunctions(elsbox.xcenter[idxbox["T"]],
                                        elsbox.ycenter[idxbox["T"]],
                                        g, rho, lambda, mu)  
    Ttx = zeros(nside)
    Tty = zeros(nside)
    for i in 1:length(Sg[:, 1])
        nvec = [elsbox.xnormal[idxbox["T"][i]] ; elsbox.ynormal[idxbox["T"][1]]]
        temp = [Sg[i, 1] Sg[i, 3] ; Sg[i, 3] Sg[i, 2]] * nvec
        Ttx[i] = temp[1]
        Tty[i] = temp[2]
    end
    bcsboxgravity[4*nside+1:2:6*nside] = Ttx
    bcsboxgravity[4*nside+2:2:6*nside] = Tty

    ###
    ### Solve the BEM problem
    ###
    bcsboxgravity *= -1
    Ueffboxparticular = THbox \ bcsboxgravity  

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
    Ueffcombined = THbox \ bcsboxgravity # FIXME: Equivalent to addind the solutions?   
    Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox,
                                   [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)

    Ugravityandfault = Ucomp .+ Uint .+ UF
    Sgravityandfault = Scomp .+ Sint .+ SF
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ugravityandfault, Sgravityandfault, "(Gravity + fault)")
    # Difference between gravity alone and gravity + fault
    # plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
    #            Ugravityonly.-Ugravityandfault, Sgravityonly.-Sgravityandfault, "(Gravity only) - (Gravity + fault)")
    # plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
    #            UF, SF, "Fullspace fault")
end
faultboxgravity()
