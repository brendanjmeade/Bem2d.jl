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
    npts = 20
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
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ubox, Sbox, "Dislocation only in box")


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
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               U, S, "Gravity on box edges only, no dislocation")


    ### This is all super-strange with a slip BC on the dislocation
    
    #
    # Box BEM problem (dislocation and gravity)
    #
    Ug, Sg = gravityparticularfunctions(elsbox.xcenter[1:1:elsbox.endidx],
                                        elsbox.ycenter[1:1:elsbox.endidx], g, rho, lambda, mu)    
    bcsboxgravity = zeros(2 * elsbox.endidx)
    bcsboxgravity[1:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 1]
    bcsboxgravity[2:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 2]
    bcsboxgravity *= -1
    bcsboxgravity[1:2] .= 0.0 # Eliminate gravity here
    Ueffboxparticular = inv(THmat) * bcsboxgravity    
    Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, collect(1:1:elsbox.endidx),
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)
    Uint, Sint = gravityparticularfunctions(xgrid, ygrid, g, rho, lambda, mu)
    Udg = @. Ucomp + Uint
    Sdg = @. Scomp + Sint
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Udg, Sdg, "Gravity on box edges and dislocation")
    @infiltrate
    
    #
    # Box BEM gravity with a dislocation
    #
    bcsboxgravity = zeros(2 * elsbox.endidx)
    bcsboxgravity[1:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 1]
    bcsboxgravity[2:2:2*idxbox["B"][end]] = Ug[1:1:idxbox["B"][end], 2]
    bcsboxgravity *= -1
    # bcsboxgravity[1:2] .+= 0.5
    bcsboxgravity[1:2] .= 0.5 # Add fault slip and eliminate gravity here
    @show cond(THmat)
    Ueffboxparticular = inv(THmat) * bcsboxgravity    
    Ucomp, Scomp = constdispstress(slip2dispstress, xgrid, ygrid, elsbox, collect(1:1:elsbox.endidx),
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)
    Uint, Sint = gravityparticularfunctions(xgrid, ygrid, g, rho, lambda, mu)
    Udg2 = @. Ucomp + Uint
    Sdg2 = @. Scomp + Sint
    plotfields(elsbox, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Udg2, Sdg2, "Gravity on box edges and dislocation")

    @infiltrate
end
dislocationinabox()
