using Revise
using DifferentialEquations
using JLD2
using Dates
using Bem2d

function calc_dvθ(vθ, p, t)
    ∂t, els, η, dc, blockvxglobal, blockvyglobal = p
    vxglobal = vθ[1:3:end]
    vyglobal = vθ[2:3:end]
    θ = vθ[3:3:end]
    vx, vy = multmatvec(els.rotmat[1:els.endidx, :, :], vxglobal, vyglobal)
    dt = ∂t * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]
    dtx, dty = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])
    dθ, dvx, dvy = zeros(els.endidx), zeros(els.endidx), zeros(els.endidx)
    for i in 1:els.endidx
        vx = abs.(vx)
        θ = abs.(θ)
        dθ[i] = -vx[i] * θ[i] / dc * log(vx[i] * θ[i] / dc) # slip law
        # dθ[i] = 1 - θ[i] * vx[i] / dc # Aging law
        dvx[i] = 1 / (η / els.σn[i] + els.a[i] / vx[i]) * (dtx[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
        dvy[i] = 0 # fault perpendicular creep
        # dvperp[i] = dvpara[i] # fault perpendicular velocity proportional to fault parallel velocity
    end
    dvxglobal, dvyglobal = multmatvec(els.rotmatinv[1:els.endidx, :, :], dvx, dvy)
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dvxglobal
    dvθ[2:3:end] = dvyglobal
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_qdmakie()
    # Fun things to play with
    nsteps = 1000
    amplitude = 1.0
    nfault = 100

    # Constants
    outfilename = string(now()) * ".jld2"
    siay = 365.25 * 24 * 60 * 60
    tspan = (0, siay * 1000)
    abstol = 1e-4
    reltol = 1e-4
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
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
        els.σn[els.endidx + i] = 50e6
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    srcidx = findall(x->x == "fault", els.name)
    obsidx = srcidx
    @time ∂u, ∂σ, ∂t = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)

    # Set initial conditions
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)
    p = (∂t, els, η, dc, blockvelx, blockvely)
    prob = ODEProblem(calc_dvθ, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol, progress = true)

    println("Step integrating")
    @time for i in 1:nsteps
        step!(integrator)
        println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
    end

    @time @save outfilename integrator.sol
    println("Wrote integration results to:")
    println(outfilename)

end
ex_qdmakie()