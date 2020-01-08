using Revise
using DifferentialEquations
using JLD2
using Dates
using Infiltrator
using PyPlot
using Bem2d

function derivsconst!(dudt, u, p, t)
    intidx, nintidx, bemsliptotractotal, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt = p
    @views multmatvec!(vx, vy, els.rotmat[intidx, :, :], u[1:3:end], u[2:3:end])
    @views dtracglobaldt = bemsliptotractotal * [blockvxglobal .- u[1:3:end] blockvyglobal .- u[2:3:end]]'[:]
    @views multmatvec!(dtracxglobaldt, dtracyglobaldt, els.rotmat[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
    for i in 1:nintidx
        @views dthetadt[i] = thetalaw(abs(vx[i]), abs(u[3:3:end][i]), dc)
        @views dvxdt[i] = 1 / (eta / els.normalstress[intidx[i]] + els.a[intidx[i]] / abs(vx[i])) * (dtracxglobaldt[i] / els.normalstress[intidx[i]] - els.b[intidx[i]] * dthetadt[i] / u[3:3:end][i])
        dvydt[i] = 0
    end
    @views multmatvec!(dudt[1:3:end], dudt[2:3:end], els.rotmatinv[intidx, :, :], dvxdt, dvydt)
    dudt[3:3:end] = dthetadt
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

function ex_qdstepandsaveplanar()
    # Constants
    norml = sqrt(2.0) / 2.0
    nsteps = 5000
    nfreesurf = 100
    nfault = 200
    printstep = 100
    amplitude = 1.0
    siay = 365.25 * 24 * 60 * 60
    tspan = (0, siay * 1000)
    abstol = 1e-4
    reltol = 1e-4
    mu = 3e10
    nu = 0.25
    rho = 2700.0
    eta = mu / (2.0 * sqrt(mu / rho))
    dc = 0.05
    blockvelx = 1e-9 * norml
    blockvely = 1e-9 * norml

    # Create fault elements
    els = Elements(Int(1e5))
    faultwidth = 10000
    
    # 45 degree fault
    x1, y1, x2, y2 = discretizedline(-5e3*norml, -5e3*norml, 5.0e3*norml, 4.0e3*norml, nfault)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.a[els.endidx + i] = 0.015
        els.b[els.endidx + i] = 0.020
        els.normalstress[els.endidx + i] = 50e7
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Create convience tools
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    
    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    @time _, _, partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    bemsliptotractotal = partialsconst["trac"]["fault"]["fault"]
    
    # CS elements - Euler style stress integration
    nnodes = 1 * nfault
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 1e-3 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)
    intidx = idx["fault"] # indices of elements to integrate
    nintidx = length(intidx)
    dthetadt = zeros(nintidx)
    dvxdt = zeros(nintidx)
    dvydt = zeros(nintidx)
    vx = zeros(nintidx)
    vy = zeros(nintidx)
    dtracxglobaldt = zeros(nintidx)
    dtracyglobaldt = zeros(nintidx)
    dtracglobaldt = zeros(2 * nintidx)
    p = (intidx, nintidx, bemsliptotractotal, els, eta, thetaaginglaw, dc, blockvelx, blockvely, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt)
    prob = ODEProblem(derivsconst!, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    @time for i in 1:nsteps
        step!(integrator)
        if mod(i, printstep) == 0
            println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
        end    
    end
    plotqdtimeseries(integrator.sol, 3, nfault)

    # Create new folder and save timestamped .jld2 file to that folder
    basename = string(now())
    outfoldername = "~/Desktop/data/qdvis/" * basename
    outfilename =  outfoldername * ".jld2"
    mkpath(outfoldername)
    @time @save outfilename integrator.sol els mu nu
    println("Wrote integration results to:")
    println(outfilename)
end
ex_qdstepandsaveplanar()
