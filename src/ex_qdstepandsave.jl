using Revise
using DifferentialEquations
using JLD2
using Dates
using Infiltrator
using Bem2d

function derivs2(u, p, t)
    partials, els, eta, dc, blockvxglobal, blockvyglobal = p
    vxglobal = u[1:5:end]
    vyglobal = u[2:5:end]
    theta = u[3:5:end]
    tx = u[4:5:end]
    ty = u[5:5:end]
    vx, vy = multmatvec(els.rotmat[1:els.endidx, :, :], vxglobal, vyglobal)
    dt =  partials["trac"]["fault"]["fault"] * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]
    dtx, dty = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])
    dtheta, dvx, dvy = zeros(els.endidx), zeros(els.endidx), zeros(els.endidx)
    vx = @. abs(vx)
    theta = @. abs(theta)

    for i in 1:els.endidx
        dtheta[i] = 1 - theta[i] * vx[i] / dc
        dvx[i] = 1 / (eta / els.normalstress[i] + els.a[i] / vx[i]) * (dtx[i] / els.normalstress[i] - els.b[i] * dtheta[i] / theta[i])
        dvy[i] = 0
    end

    dvxglobal, dvyglobal = multmatvec(els.rotmatinv[1:els.endidx, :, :], dvx, dvy)
    dudt = zeros(5 * els.endidx)
    dudt[1:5:end] = dvxglobal
    dudt[2:5:end] = dvyglobal
    dudt[3:5:end] = dtheta
    dudt[4:5:end] = dtx
    dudt[5:5:end] = dty
    return dudt
end

function derivs(u, p, t)
    partials, els, eta, dc, blockvxglobal, blockvyglobal = p
    vxglobal = u[1:3:end]
    vyglobal = u[2:3:end]
    theta = u[3:3:end]
    vx, vy = multmatvec(els.rotmat[1:els.endidx, :, :], vxglobal, vyglobal)
    dt =  partials["trac"]["fault"]["fault"] * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]
    dtx, dty = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])
    dtheta, dvx, dvy = zeros(els.endidx), zeros(els.endidx), zeros(els.endidx)
    vx = @. abs(vx)
    theta = @. abs(theta)

    for i in 1:els.endidx
        dtheta[i] = 1 - theta[i] * vx[i] / dc
        dvx[i] = 1 / (eta / els.normalstress[i] + els.a[i] / vx[i]) * (dtx[i] / els.normalstress[i] - els.b[i] * dtheta[i] / theta[i])
        dvy[i] = 0
    end

    dvxglobal, dvyglobal = multmatvec(els.rotmatinv[1:els.endidx, :, :], dvx, dvy)
    dudt = zeros(3 * els.endidx)
    dudt[1:3:end] = dvxglobal
    dudt[2:3:end] = dvyglobal
    dudt[3:3:end] = dtheta
    return dudt
end

function derivsquad(u, p, t)
    partials, els, eta, dc, blockvxglobal, blockvyglobal = p
    vxglobal = u[1:3:end]
    vyglobal = u[2:3:end]
    theta = u[3:3:end]
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
    @time _, _, partialsquad["trac"]["fault"]["fault"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)

    #
    # CS elements - Euler style stress integration
    #
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)
    p = (partialsconst, els, eta, dc, blockvelx, blockvely)
    prob = ODEProblem(derivs, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    @time for i in 1:nsteps
        step!(integrator)
        # println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
    end
    @show integrator.sol.t
    # plotqdtimeseries(integrator.sol, 3, nfault)

    
    #
    # 3QN elements - Euler style stress integration
    #
    nnodes = 3 * nfault
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)
    p = (partialsquad, els, eta, dc, blockvelx, blockvely)
    prob = ODEProblem(derivsquad, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    @time for i in 1:nsteps
        step!(integrator)
        # println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
    end
    @show integrator.sol.t
    # plotqdtimeseries(integrator.sol, 3, nfault)

    
    
    #
    # Stresses fully incorporated in state vector
    #
    # ics = zeros(5 * nnodes)
    # ics[1:5:end] = 1e-3 * blockvelx * ones(nnodes)
    # ics[2:5:end] = 0.0 * blockvely * ones(nnodes)
    # ics[3:5:end] = 1e8 * ones(nnodes)
    # ics[4:5:end] = 0.0 * ones(nnodes)
    # ics[5:5:end] = 0.0 * ones(nnodes)
    # p = (partialsconst, els, eta, dc, blockvelx, blockvely)
    # prob = ODEProblem(derivs2, ics, tspan, p)
    # integrator2 = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    # @time for i in 1:nsteps
    #     step!(integrator2)
    #     # println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
    # end
    # plotqdtimeseries(integrator2.sol, 5, nfault)

    # @time @save outfilename integrator.sol
    # println("Wrote integration results to:")
    # println(outfilename)
end
ex_qdstepandsave()
