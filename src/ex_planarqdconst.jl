using Revise
using DifferentialEquations
# using PyCall
# using PyPlot
# using Colors
# using ColorSchemes
using JLD2
using Dates
using UUIDs
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
        vx = abs.(vx) # Chris R. suggestion
        θ = abs.(θ) # Chris R. suggestion
        dθ[i] = -vx[i] * θ[i] / dc * log(vx[i] * θ[i] / dc) # slip law
        # dθ[i] = 1 - θ[i] * vx[i] / dc # Aging law
        dvx[i] = 1 / (η / els.σn[i] + els.a[i] / vx[i]) * (dtx[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
        dvy[i] = 0 # fault perpendicular creep
        # dvperp[i] = dvpara[i] # fault perpendicular velocity proportional to fault parallel velocity
    end
    dvxglobal, dvyglobal = Bem2d.multmatvec(els.rotmatinv[1:els.endidx, :, :], dvx, dvy)
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dvxglobal
    dvθ[2:3:end] = dvyglobal
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_planarqdconst()
    # Constants and model parameters
    outfilename = string(Dates.now()) * '_' * string(UUIDs.uuid4()) * ".jld2"
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
    els = Bem2d.Elements(Int(1e5))
    nfault = 30
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = Bem2d.discretizedline(-faultwidth, 0, faultwidth, 0, nfault)

    # Modify y1, and y2 for a sinusoidal fault
    amplitude = 1000.0
    y1 = amplitude * @.sin(2 * π * x1 / faultwidth)
    y2 = amplitude * @.sin(2 * π * x2 / faultwidth)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.a[els.endidx + 1:els.endidx + nfault] .= 0.015
    els.b[els.endidx + 1:els.endidx + nfault] .= 0.020
    els.σn[els.endidx + 1:els.endidx + nfault] .= 50e6
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    Bem2d.standardize_elements!(els)

    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    srcidx = findall(x->x == "fault", els.name)
    obsidx = srcidx
    @time ∂u, ∂σ, ∂t = Bem2d.∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)

    # Set initial conditions
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)

    # (Bulk) Time integrate elastic model
    p = (∂t, els, η, dc, blockvelx, blockvely)
    prob = DifferentialEquations.ODEProblem(calc_dvθ, ics, tspan, p)
    @time sol = solve(prob, DifferentialEquations.DP5(), abstol = abstol, reltol = reltol, progress = true)
    @time Bem2d.plottimeseries(sol)
    @time @save outfilename sol
    println("Wrote integration results to:")
    println(outfilename)
end
ex_planarqdconst()
