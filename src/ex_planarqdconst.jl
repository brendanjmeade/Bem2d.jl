using Revise
using DifferentialEquations
using PyCall
using PyPlot
using Bem2d

# Derivatives to feed to ODE integrator
function calc_dvθ(vθ, p, t)
    display(t)
    ∂t, els, η, dc = p
    θ = vθ[3:3:end]
    v = @.sqrt(vθ[1:3:end].^2 + vθ[2:3:end].^2) # TODO: Flat fault only
    dt = ∂t * [vθ[1:3:end] ; vθ[2:3:end]]
    dτ = dt[1:2:end]

    # State evolution: aging law
    dθ = -v .* θ ./ dc .* @.log(v .* θ ./ dc)

    # Acceleration evolution
    dv = 1 ./ (η ./ els.σn[1:els.endidx] .+ els.a[1:els.endidx] ./ v) .*
        (dτ ./ els.σn[1:els.endidx] .- els.b[1:els.endidx] .* dθ ./ θ)
    
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] .= dv # TODO: Flat fault only
    dvθ[2:3:end] .= 0 # TODO: Flat fault only
    dvθ[3:3:end] .= dθ
    return dvθ
end

function ex_planarqdconst()
    # Constants and model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 1000.00)
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    dc = 0.2
    blockvelx = 1e-9
    blockvely = 0.0

    # Create fault elements
    els = Elements()
    nfault = 50
    nnodes = 1 * nfault
    L = 10000
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nfault)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.a[els.endidx + 1:els.endidx + nfault] .= 0.015
    els.b[els.endidx + 1:els.endidx + nfault] .= 0.020
    els.σn[els.endidx + 1:els.endidx + nfault] .= 30e6
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    srcidx = findall(x->x == "fault", els.name)
    obsidx = srcidx
    _, _, ∂t = ∂constslip(els, srcidx, obsidx, μ, ν)

    # Set initial conditions
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 1e-3 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)

    # Time integrate
    p = (∂t, els, η, dc)
    prob = ODEProblem(calc_dvθ, ics, tspan, p)
    sol = solve(prob, abstol = 1e-4, reltol = 1e-4, StepsizeLimiter = false)
    
    close("all")
    plot(sol.t, "-b")
    show()
end
ex_planarqdconst()
