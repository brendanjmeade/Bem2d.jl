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
    xgrid, ygrid = obsgrid(-30e3, -20e3, 30e3, -1, npts)

    # Element geometries and data structures for half space approximation
    elshalfspace = Elements(Int(1e5))
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(elshalfspace, x1, y1, x2, y2, "fault")
    nsurf = 100
    x1, y1, x2, y2 = discretizedline(-100e3, 0, 100e3, 0, nsurf) # Free surface
    addelsez!(elshalfspace, x1, y1, x2, y2, "surf")
    idxhalfspace = getidxdict(elshalfspace)

    # Element geometries and data structures for the box cases
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
    # Halfspace BEM problem
    #
    _, H_surf_fault = PUTC(slip2dispstress, elshalfspace,
                           idxhalfspace["surf"], idxhalfspace["fault"], mu, nu)
    _, H_surf_surf = PUTC(slip2dispstress, elshalfspace,
                          idxhalfspace["surf"], idxhalfspace["surf"], mu, nu)
    faultslip = sqrt(2) / 2 * ones(2 * nfault)
    Ueff = inv(H_surf_surf) * (H_surf_fault * faultslip)

    # Forward evaluation
    Ufault, Sfault = constdispstress(slip2dispstress, xgrid, ygrid, elshalfspace,
                                     idxhalfspace["fault"], ones(nfault), ones(nfault),
                                     mu, nu)
    Usurf, Ssurf = constdispstress(slip2dispstress, xgrid, ygrid, elshalfspace,
                                   idxhalfspace["surf"], Ueff[1:2:end], Ueff[2:2:end],
                                   mu, nu)
    Utotal = @. Ufault - Usurf
    Stotal = @. Sfault - Ssurf
    plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Utotal, Stotal, "Total")

    # Alternative BEM solution
    T_fault_all, _ = PUTC(slip2dispstress, elshalfspace, idxhalfspace["fault"],
                          collect(1:1:elshalfspace.endidx), mu, nu)
    _, H_surf_all = PUTC(slip2dispstress, elshalfspace, idxhalfspace["surf"],
                         collect(1:1:elshalfspace.endidx), mu, nu)
    bcs = zeros(2 * elshalfspace.endidx)
    bcs[1:2*nfault] .= 0.5
    Ueffalt = inv([T_fault_all; alpha .* H_surf_all]) * bcs
    
    # Alternative forward evaluation
    Ualt, Salt = constdispstress(slip2dispstress, xgrid, ygrid, elshalfspace,
        collect(1:1:elshalfspace.endidx), Ueffalt[1:2:end], Ueffalt[2:2:end], mu, nu)
    plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ualt, Salt, "Alternative")

    # Plot difference
    Udiff = @. Ualt - Utotal
    Sdiff = @. Salt - Stotal
    plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Udiff, Sdiff, "Alternative - Total")

    # Okada comparison if on Linux.
    # I can only PyCall to see okada_wrapper there :(
    if Sys.islinux()
        ow = pyimport("okada_wrapper") # from okada_wrapper import dc3dwrapper
        Uokada = zeros(length(xgrid), 2)
        Sokada = zeros(length(xgrid), 3)
        for i in 1:length(xgrid)
            _, u, gradu = ow.dc3dwrapper((2.0/3.0),
                                         [0.0, xgrid[i]+5e3, ygrid[i]],
                                         5e3,
                                         45,
                                         [-1000000, 1000000],
                                         [-5e3*sqrt(2), 5e3*sqrt(2)],
                                         [0.0, 1.0, 0.0])
            strain = @. 0.5 * (gradu' + gradu)
            stress = mu*LinearAlgebra.I(3)*tr(strain) + 2.0*mu*strain
            Uokada[i, 1] = u[2]
            Uokada[i, 2] = u[3]
            Sokada[i, 1] = stress[2, 2]
            Sokada[i, 2] = stress[3, 3]
            Sokada[i, 3] = stress[2, 3]
        end
        plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
                   Uokada, Sokada, "Okada")

        # Let's look at some Residuals
        Utotalres = @. Uokada - Utotal
        Stotalres = @. Sokada - Stotal
        Ualtres = @. Uokada - Ualt
        Saltres = @. Sokada - Salt
        plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Utotalres, Stotalres, "Okada - Total")
        plotfields(elshalfspace, reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
               Ualtres, Saltres, "Okada - Alternative")
    end

    #
    # Box BEM problem
    #


end
dislocationinabox()
