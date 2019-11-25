using Revise
using DifferentialEquations
using JLD2
using Dates
using Infiltrator
using Bem2d

function derivsconst(u, p, t)
    partials, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal = p
    vxglobal = @. abs(u[1:3:end])
    vyglobal = u[2:3:end]
    theta = @. abs(u[3:3:end])
    dtracglobaldt =  partials["trac"]["fault"]["fault"] * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]

    vx, vy = multmatvec(els.rotmat[1:els.endidx, :, :], vxglobal, vyglobal)
    dtracxglobaldt, dtracyglobaldt = multmatvec(els.rotmat[1:els.endidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])

    dthetadt = zeros(els.endidx)
    dvxdt = zeros(els.endidx)
    dvydt = zeros(els.endidx)    
    for i in 1:els.endidx
        # dthetadt[i] = 1 - theta[i] * vx[i] / dc
        dthetadt[i] = thetalaw(vx[i], theta[i], dc)
        dvxdt[i] = 1 / (eta / els.normalstress[i] + els.a[i] / vx[i]) * (dtracxglobaldt[i] / els.normalstress[i] - els.b[i] * dthetadt[i] / theta[i])
        dvydt[i] = 0
    end

    dvxglobaldt, dvyglobaldt = multmatvec(els.rotmatinv[1:els.endidx, :, :], dvxdt, dvydt)
    dudt = zeros(3 * els.endidx)
    dudt[1:3:end] = dvxglobaldt
    dudt[2:3:end] = dvyglobaldt
    dudt[3:3:end] = dthetadt
    return dudt
end

function derivsquad(u, p, t)
    partials, els, eta, dc, blockvxglobal, blockvyglobal = p
    vxglobal = u[1:3:end]
    vyglobal = u[2:3:end]
    theta = u[3:3:end]
    @infiltrate
    println("here - 1")
    vx, vy = multmatvec(repeat(els.rotmat[1:els.endidx, :, :], 3, 1, 1), vxglobal, vyglobal)
    println("here - 2")
    dt =  partials["trac"]["fault"]["fault"] * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]
    println("here - 3")
    dtx, dty = multmatvec(repeat(els.rotmat[1:els.endidx, :, :], 3, 1, 1), dt[1:2:end], dt[2:2:end])
    println("here - 4")
    dtheta, dvx, dvy = zeros(3*els.endidx), zeros(3*els.endidx), zeros(3*els.endidx)
    vx = @. abs(vx)
    theta = @. abs(theta)

    for i in 1:3*els.endidx
        dtheta[i] = 1 - theta[i] * vx[i] / dc
        dvx[i] = 1 / (eta / els.normalstress[i] + els.a[i] / vx[i]) * (dtx[i] / els.normalstress[i] - els.b[i] * dtheta[i] / theta[i])
        dvy[i] = 0
    end

    dvxglobal, dvyglobal = multmatvec(repeat(els.rotmatinv[1:els.endidx, :, :], 3, 1, 1), dvx, dvy)
    dudt = zeros(9 * els.endidx)
    dudt[1:3:end] = dvxglobal
    dudt[2:3:end] = dvyglobal
    dudt[3:3:end] = dtheta
    return dudt
end


function ex_qdstepandsave()
    # Constants
    nsteps = 500
    printstep = 100
    amplitude = 1.0
    nfault = 100
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
    els = Elements(Int(1e5))
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = discretizedline(-faultwidth, 0, faultwidth, 0, nfault)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.a[els.endidx + i] = 0.015
        els.b[els.endidx + i] = 0.020
        els.normalstress[els.endidx + i] = 50e6
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Create convience tools
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    partialsquad = initpartials(els)
    
    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    @time _, _, partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    # @time _, _, partialsquad["trac"]["fault"]["fault"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)

    #
    # CS elements - Euler style stress integration
    #
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)
    p = (partialsconst, els, eta, thetaaginglaw, dc, blockvelx, blockvely)
    prob = ODEProblem(derivsconst, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    @time for i in 1:nsteps
        step!(integrator)
        if mod(i, printstep) == 0
            println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
        end    
    end
    # @show integrator.sol.t
    plotqdtimeseries(integrator.sol, 3, nfault)

    
    #
    # 3QN elements - Euler style stress integration
    #
    # nnodes = 3 * nfault
    # ics = zeros(3 * nnodes)
    # ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    # ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    # ics[3:3:end] = 1e8 * ones(nnodes)
    # p = (partialsquad, els, eta, dc, blockvelx, blockvely)
    # prob = ODEProblem(derivsquad, ics, tspan, p)
    # integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    # @time for i in 1:nsteps
    #     step!(integrator)
    #     # println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
    # end
    # @show integrator.sol.t
    # plotqdtimeseries(integrator.sol, 3, nfault)

    # @time @save outfilename integrator.sol
    # println("Wrote integration results to:")
    # println(outfilename)
end
ex_qdstepandsave()
