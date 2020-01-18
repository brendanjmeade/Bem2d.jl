using Revise
using DifferentialEquations
using JLD2
using Dates
using Infiltrator
using PyPlot
using Bem2d

# function derivsconst(u, p, t)
#     intidx, partials, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal = p
#     nintidx = length(intidx)
#     vxglobal = @. abs(u[1:3:end])
#     vyglobal = u[2:3:end]
#     theta = @. abs(u[3:3:end])
#     dtracglobaldt =  partials["trac"]["fault"]["fault"] * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]
#     vx, vy = multmatvec(els.rotmatinv[intidx, :, :], vxglobal, vyglobal)
#     dtracxglobaldt, dtracyglobaldt = multmatvec(els.rotmatinv[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
    
#     dthetadt = zeros(nintidx)
#     dvxdt = zeros(nintidx)
#     dvydt = zeros(nintidx)
#     for i in 1:length(intidx)
#         dthetadt[i] = thetalaw(vx[i], theta[i], dc)
#         dvxdt[i] = 1 / (eta / els.normalstress[intidx[i]] + els.a[intidx[i]] / vx[i]) * (dtracxglobaldt[i] / els.normalstress[intidx[i]] - els.b[intidx[i]] * dthetadt[i] / theta[i])
#         dvydt[i] = 0
#     end

#     dvxglobaldt, dvyglobaldt = multmatvec(els.rotmat[intidx, :, :], dvxdt, dvydt)
#     dudt = zeros(length(u))
#     dudt[1:3:end] = dvxglobaldt
#     dudt[2:3:end] = dvyglobaldt
#     dudt[3:3:end] = dthetadt
#     return dudt
# end

function derivsconst!(dudt, u, p, t)
    intidx, nintidx, bemsliptotractotal, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt = p
    @views multmatvec!(vx, vy, els.rotmatinv[intidx, :, :], u[1:3:end], u[2:3:end])
    @views dtracglobaldt = bemsliptotractotal * [blockvxglobal .- u[1:3:end] blockvyglobal .- u[2:3:end]]'[:]
    @views multmatvec!(dtracxglobaldt, dtracyglobaldt, els.rotmatinv[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
    for i in 1:nintidx
        @views dthetadt[i] = thetalaw(abs(vx[i]), abs(u[3:3:end][i]), dc)
        @views dvxdt[i] = 1 / (eta / els.normalstress[intidx[i]] + els.a[intidx[i]] / abs(vx[i])) * (dtracxglobaldt[i] / els.normalstress[intidx[i]] - els.b[intidx[i]] * dthetadt[i] / u[3:3:end][i])
        dvydt[i] = 0
    end
    @views multmatvec!(dudt[1:3:end], dudt[2:3:end], els.rotmat[intidx, :, :], dvxdt, dvydt)
    dudt[3:3:end] = dthetadt
    return nothing
end

# function derivsconst!(dudt, u, p, t)
#     intidx, nintidx, bemsliptotractotal, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt = p
#     @views multmatvec!(vx, vy, els.rotmat[intidx, :, :], u[1:3:end], u[2:3:end])
#     @views dtracglobaldt = bemsliptotractotal * [blockvxglobal .- u[1:3:end] blockvyglobal .- u[2:3:end]]'[:]
#     @views multmatvec!(dtracxglobaldt, dtracyglobaldt, els.rotmat[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
#     for i in 1:nintidx
#         @views dthetadt[i] = thetalaw(abs(vx[i]), abs(u[3:3:end][i]), dc)
#         @views dvxdt[i] = 1 / (eta / els.normalstress[intidx[i]] + els.a[intidx[i]] / abs(vx[i])) * (dtracxglobaldt[i] / els.normalstress[intidx[i]] - els.b[intidx[i]] * dthetadt[i] / u[3:3:end][i])
#         dvydt[i] = 0
#     end
#     @views multmatvec!(dudt[1:3:end], dudt[2:3:end], els.rotmatinv[intidx, :, :], dvxdt, dvydt)
#     dudt[3:3:end] = dthetadt
#     return nothing
# end

function derivsquad!(dudt, u, p, t)
    intidx, nintidx, partials, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt = p
    @views multmatvecquad!(vx, vy, els.rotmatinv[intidx, :, :], u[1:3:end], u[2:3:end])
    @views dtracglobaldt = partials["trac"]["fault"]["fault"] * [blockvxglobal .- u[1:3:end] blockvyglobal .- u[2:3:end]]'[:]
    @views multmatvecquad!(dtracxglobaldt, dtracyglobaldt, els.rotmatinv[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
    for i in 1:3*nintidx
        elidx = Int64(floor((i - 1) / 3) + 1) # Change w/ every 3rd node
        @views dthetadt[i] = thetalaw(abs(vx[i]), abs(u[3:3:end][i]), dc)
        @views dvxdt[i] = 1 / (eta / els.normalstress[intidx[elidx]] + els.a[intidx[elidx]] / abs(vx[i])) * (dtracxglobaldt[i] / els.normalstress[intidx[elidx]] - els.b[intidx[elidx]] * dthetadt[i] / u[3:3:end][i])
        dvydt[i] = 0
    end
    @views multmatvecquad!(dudt[1:3:end], dudt[2:3:end], els.rotmat[intidx, :, :], dvxdt, dvydt)
    dudt[3:3:end] = dthetadt
    return nothing
end

function ex_qdstepandsave()
    # Constants
    nsteps = 5000
    nfreesurf = 100
    nfault = 200
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
    els = Elements(Int(1e5))
    faultwidth = 10000
    
    # Curved fault
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, nfault)
    y1 = 3e3 * atan.(x1 / 1e3)
    y2 = 3e3 * atan.(x2 / 1e3)
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

    # Topographic free surface
    x1, y1, x2, y2 = discretizedline(-20e3, 0, 20e3, 0, nfreesurf)
    y1 = -1e3 * atan.(x1 / 1e3)
    y2 = -1e3 * atan.(x2 / 1e3)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "freesurftopo"
        els.a[els.endidx + i] = 0.015
        els.b[els.endidx + i] = 0.020
        els.normalstress[els.endidx + i] = 0e6
    end
    standardize_elements!(els)
    

    # Create convience tools
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    
    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    @time _, _, partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    @time _, _, partialsconst["trac"]["fault"]["freesurftopo"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurftopo"], mu, nu)
    @time _, _, partialsconst["trac"]["freesurftopo"]["freesurftopo"] = partialsconstdispstress(slip2dispstress, els, idx["freesurftopo"], idx["freesurftopo"], mu, nu)
    @time _, _, partialsconst["trac"]["freesurftopo"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["freesurftopo"], idx["fault"], mu, nu)
    
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
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
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
    prob = DifferentialEquations.ODEProblem(derivsconst!, ics, tspan, p)
    integrator = DifferentialEquations.init(prob, Vern7(), abstol = abstol, reltol = reltol)
    @time for i in 1:nsteps
        DifferentialEquations.step!(integrator)
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
