using Revise
using LinearAlgebra
using DifferentialEquations
using JLD2
using Dates
using Infiltrator
using PyPlot
using Bem2d

# Have access to 500 xeon cores
# Dedicated access to:
# * One server has 3TB RAM, 8 P100 GPUs, 80 cores
# * Two servers have 1TB RAM, 32 cores, 8 V100 GPUs

function derivsconst!(dudt, u, p, t)
    bemsliptotractotal, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt = p
    # @views multmatvec!(vx, vy, els.rotmat[intidx, :, :], u[1:3:end], u[2:3:end])
    @inbounds for i in 1:length(els)
        vx[i] = els[i].rotmat[1, 1] * u[(1:3:end)[i]] + els[i].rotmat[1, 2] * u[(2:3:end)[i]]
        vy[i] = els[i].rotmat[2, 1] * u[(1:3:end)[i]] + els[i].rotmat[2, 2] * u[(2:3:end)[i]]
    end

    # THIS IS WHERE MOST TIME IS SPENT:
    mul!(dtracglobaldt, bemsliptotractotal, @views([blockvxglobal .- u[1:3:end] blockvyglobal .- u[2:3:end]]')[:])
    # @views multmatvec!(dtracxglobaldt, dtracyglobaldt, els.rotmat[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
    @inbounds for i in 1:length(els)
        dtracxglobaldt[i] = els[i].rotmat[1, 1] * dtracglobaldt[(1:2:end)[i]] + els[i].rotmat[1, 2] * dtracglobaldt[(2:2:end)[i]]
        dtracyglobaldt[i] = els[i].rotmat[2, 1] * dtracglobaldt[(1:2:end)[i]] + els[i].rotmat[2, 2] * dtracglobaldt[(2:2:end)[i]]
    end
    for i in 1:length(els)
        dthetadt = thetalaw(abs(vx[i]), abs(u[(3:3:end)[i]]), dc)
        dvxdt[i] = 1 / (eta / els[i].normalstress + els[i].a / abs(vx[i])) * (dtracxglobaldt[i] / els[i].normalstress - els[i].b * dthetadt / u[(3:3:end)[i]])
        dvydt[i] = 0
        dudt[(3:3:end)[i]] = dthetadt
    end
    # @views multmatvec!(dudt[1:3:end], dudt[2:3:end], els.rotmatinv[intidx, :, :], dvxdt, dvydt)
    @inbounds for i in 1:length(els)
        dudt[(1:3:end)[i]] = els[i].rotmatinv[1, 1] * dvxdt[i] + els[i].rotmatinv[1, 2] * dvydt[i]
        dudt[(2:3:end)[i]] = els[i].rotmatinv[2, 1] * dvxdt[i] + els[i].rotmatinv[2, 2] * dvydt[i]
    end
    return nothing
end

function derivsquad!(dudt, u, p, t)
    intidx, nintidx, partials, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt = p
    @views multmatvecquad!(vx, vy, els.rotmat[intidx, :, :], u[1:3:end], u[2:3:end])
    @views dtracglobaldt = partials["trac"]["fault"]["fault"] * [blockvxglobal .- u[1:3:end] blockvyglobal .- u[2:3:end]]'[:]
    @views multmatvecquad!(dtracxglobaldt, dtracyglobaldt, els.rotmat[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
    for i in 1:3*nintidx
        elidx = Int64(floor((i - 1) / 3) + 1) # Change w/ every 3rd node
        @views dthetadt[i] = thetalaw(abs(vx[i]), abs(u[3:3:end][i]), dc)
        @views dvxdt[i] = 1 / (eta / els.normalstress[intidx[elidx]] + els.a[intidx[elidx]] / abs(vx[i])) * (dtracxglobaldt[i] / els.normalstress[intidx[elidx]] - els.b[intidx[elidx]] * dthetadt[i] / u[3:3:end][i])
        dvydt[i] = 0
    end
    @views multmatvecquad!(dudt[1:3:end], dudt[2:3:end], els.rotmatinv[intidx, :, :], dvxdt, dvydt)
    dudt[3:3:end] = dthetadt
    return nothing
end

using Serialization

function ex_qdstepandsave()
    # Constants
    nsteps = 5000
    nfreesurf = 100
    nfault = 200 # can increase to 100k. When done single-threaded on CPU, take 3-10 days.
    printstep = 100
    amplitude = 1.0
    outfilename = string(now()) * ".jld2"
    siay = 365.25 * 24 * 60 * 60
    tspan = (0, siay * 1000)
    abstol = 1e-4
    reltol = 1e-4
    mu = 3e10
    nu = 0.25
    rho = 2700.0
    eta = mu / (2.0 * sqrt(mu / rho))
    dc = 0.05
    blockvelx = 1e-9
    blockvely = 0.0

    # Create fault elements
    faultwidth = 10000

    # Curved fault
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, nfault)
    y1 = 3e3 .* atan.(x1 ./ 1e3)
    y2 = 3e3 .* atan.(x2 ./ 1e3)
    fault = [
        Element(
            x1=x1[i],
            y1=y1[i],
            x2=x2[i],
            y2=y2[i],
            a=0.015,
            b=0.020,
            normalstress=50e7)
        for i in 1:length(x1)
    ]

    # Topographic free surface
    x1, y1, x2, y2 = discretizedline(-20e3, 0, 20e3, 0, nfreesurf)
    y1 = -1e3 .* atan.(x1 ./ 1e3)
    y2 = -1e3 .* atan.(x2 ./ 1e3)
    freesurftopo = [
        Element(
            x1 = x1[i],
            y1 = y1[i],
            x2 = x2[i],
            y2 = y2[i],
            a = 0.015,
            b = 0.020,
            normalstress = 0e60)
        for i in 1:length(x1)
    ]

    # partials (wierd name for partial derviatives/setup)
    # const (as opposed to a quadratic)
    partialsconst = initpartials("fault", "freesurftopo")

    # Calculate slip to traction partials on the fault
    # These are expensive.
    println("Calculating velocity to traction matrix")
    @time _, _, partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, fault, fault, mu, nu)
    @time _, _, partialsconst["trac"]["fault"]["freesurftopo"] = partialsconstdispstress(slip2dispstress, fault, freesurftopo, mu, nu)
    @time _, _, partialsconst["trac"]["freesurftopo"]["freesurftopo"] = partialsconstdispstress(slip2dispstress, freesurftopo, freesurftopo, mu, nu)
    @time _, _, partialsconst["trac"]["freesurftopo"]["fault"] = partialsconstdispstress(slip2dispstress, freesurftopo, fault, mu, nu)

    #
    # Just trying out some new notation...need to think about it.
    # parC["T"]["fault"]["fault"]
    # parQ["S"]["fault"]["fault"]
    #

    #
    # Solve the BEM problem once so that it can be passed to the solver
    # and we only have to do the matrix-vector multiply in solver
    # Note that we need both the fault and surface contribution to stress
    # on the fault.
    #

    # This matrix allows us to go from fault slip to surface displacements
    # We then need to go from surface displacements to tractions on the fault
    # Do we need tractions on the surface too to go to fault tractions?  No because it's a free surface.
    @time faultsliptosurfacedisp = inv(partialsconst["trac"]["freesurftopo"]["freesurftopo"]) * partialsconst["trac"]["fault"]["freesurftopo"]
    faultsliptosurfacedisptofaulttraction = partialsconst["trac"]["freesurftopo"]["fault"] * faultsliptosurfacedisp
    bemsliptotractotal = partialsconst["trac"]["fault"]["fault"] - faultsliptosurfacedisptofaulttraction

    # CS elements - Euler style stress integration
    nnodes = 1 * nfault
    ics = zeros(3 * nnodes) # ics: initial conditions
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)
    intels = fault # elements to integrate
    # scratch buffers to use during the ODE solve
    nintels = length(intels)
    dthetadt = zeros(nintels)
    dvxdt = zeros(nintels)
    dvydt = zeros(nintels)
    vx = zeros(nintels)
    vy = zeros(nintels)
    dtracxglobaldt = zeros(nintels)
    dtracyglobaldt = zeros(nintels)
    dtracglobaldt = zeros(2 * nintels)
    p = (bemsliptotractotal, fault, eta, thetaaginglaw, dc, blockvelx, blockvely, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt)
    prob = ODEProblem(derivsconst!, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    @time for i in 1:nsteps
        # Do this iteratively for monitoring — some runs take weeks.
        # Use Makie graphics here sometimes to see what's happening
        step!(integrator)
        if mod(i, printstep) == 0
            println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
        end
    end
    plotqdtimeseries(integrator.sol, 3, nfault)

    # 3QN elements - Euler style stress integration
    # nnodes = 3 * nfault
    # ics = zeros(3 * nnodes)
    # ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    # ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    # ics[3:3:end] = 1e8 * ones(nnodes)
    # intidx = collect(1:1:els.endidx) # indices of elements to integrate
    # nintidx = length(intidx)
    # dthetadt = zeros(3 * nintidx)
    # dvxdt = zeros(3 * nintidx)
    # dvydt = zeros(3 * nintidx)
    # vx = zeros(3 * nintidx)
    # vy = zeros(3 * nintidx)
    # dtracxglobaldt = zeros(3 * nintidx)
    # dtracyglobaldt = zeros(3 * nintidx)
    # dtracglobaldt = zeros(2 * 3 * nintidx)
    # p = (intidx, nintidx, partialsquad, els, eta, thetasliplaw, dc, blockvelx, blockvely, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt)
    # prob = ODEProblem(derivsquad!, ics, tspan, p)
    # integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    # @time for i in 1:nsteps
    #     step!(integrator)
    #     if mod(i, printstep) == 0
    #         println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
    #     end
    # end
    # plotqdtimeseriesquad(integrator.sol, 3, nfault)

    @time @save outfilename integrator.sol els mu nu
    println("Wrote integration results to:")
    println(outfilename)
end
ex_qdstepandsave()
