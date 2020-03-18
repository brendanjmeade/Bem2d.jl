using Revise
using DifferentialEquations
using JLD2
using Dates
using Infiltrator
using PyPlot
using Bem2d

function stressfunction(slip, mu, x, y, d1, d2)
    stressd1 = @. slip*mu/(2*pi) * ((y+d1)/(x^2+(y+d1)^2) - (y+d1)/(x^2+(y-d1)^2))
    stressd2 = @. slip*mu/(2*pi) * ((y+d2)/(x^2+(y+d2)^2) - (y+d2)/(x^2+(y-d2)^2))
    stress = stressd1 .- stressd2
    return stress
end

function interactionmatrix(nels, dtop, dbot, mu)
    mat = zeros(nels, nels) # to store results
    dvec = collect(LinRange(dtop, dbot, nels+1))
    dtopvec = dvec[1:1:end-1]
    dbotvec = dvec[2:1:end]
    ymid = @. -(dtopvec + dbotvec) / 2 
    
    # Loop over all possible combinations
    for i = 1:nels
        for j = 1:nels
            mat[i, j] = stressfunction(1, mu, 0, ymid[i], dtopvec[j], dbotvec[j])
        end
    end
    return mat
end

function derivsconst!(dudt, u, p, t)
    nels, U2Tmat, eta, a, b, sigma, thetalaw, dc, blockvel, dthetadt, dvdt, v, dTdt = p
    @views dTdt = U2Tmat * (blockvel .- u[1:2:end])
    for i in 1:nels
        @views dthetadt[i] = thetalaw(abs(v[i]), abs(u[2:2:end][i]), dc)
        @views dvdt[i] = 1 / (eta / sigma[i] + a[i] / abs(v[i])) * (dTdt[i] / sigma[i] - b[i] * dthetadt[i] / u[1:2:end][i])
    end
    dudt[1:2:end] = dvdt
    dudt[2:2:end] = dthetadt
    return nothing
end

function strikeslipstress()
    close("all")

    # Constants to calculate interaction matrix
    nels = 100
    mu = 3e10
    mindepth = 0e3 # meters
    maxdepth = 20e3 # meters

    #! Constants for QD
    nsteps = 500
    printstep = 100
    outfilename = string(now()) * ".jld2"
    siay = 365.25 * 24 * 60 * 60
    tspan = (0, siay * 1000)
    abstol = 1e-4
    reltol = 1e-4
    nu = 0.25
    rho = 2700.0
    eta = mu / (2.0 * sqrt(mu / rho))
    dc = 0.05
    blockvel = 1e-9
    a = 0.015 * ones(nels)
    b = 0.020 * ones(nels)
    sigma = 50e6 * ones(nels)

    #! Calculate full element to element interaction interaction
    @time U2Tmat = interactionmatrix(nels, mindepth, maxdepth, mu)
    matshow(log10.(abs.(U2Tmat))) # Plot interaction matrix

    #! Set initial conditions and package parameters for pass to integrator
    ics = zeros(2 * nels)
    ics[1:2:end] = 1e-3 * blockvel * ones(nels)
    ics[2:2:end] = 1e8 * ones(nels)
    dthetadt = zeros(nels)
    dvdt = zeros(nels)
    v = zeros(nels)
    dTdt = zeros(nels)

    #! Set up integrator
    p = (nels, U2Tmat, eta, a, b, sigma, thetaaginglaw, dc, blockvel, dthetadt, dvdt, v, dTdt)
    prob = ODEProblem(derivsconst!, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)

    #! Time integrate
    @time for i in 1:nsteps
        step!(integrator)
        if mod(i, printstep) == 0
            println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
        end    
    end

    #! Plot and save
    # plotqdtimeseries(integrator.sol, 2, nfault)    
    # @time @save outfilename integrator.sol mu nu
    # println("Wrote integration results to:")
    # println(outfilename)
end
strikeslipstress()
