using Revise
using DifferentialEquations
using PyCall
using PyPlot
using Bem2d

# Derivatives to feed to ODE integrator
function calc_dvθ(vθ, p, t)
    display(t / (365.25 * 24 * 60 * 60))
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
    dvθ[1:3:end] = dv # TODO: Flat fault only
    dvθ[2:3:end] .= 0 # TODO: Flat fault only
    dvθ[3:3:end] = dθ
    return dvθ
end


# 1-D Derivatives to feed to ODE integrator
function calc_dvθ1d(vθ, p, t)
    display(t / (365.25 * 24 * 60 * 60))
    ∂t, els, η, dc = p
    θ = vθ[3:3:end]
    v = vθ[1:3:end]

    # Hacky "1d" version
    a = 0.015
    b = 0.02
    σn = 30e6
    μ = 3e10
    Vp = 1e-9
    dτ = μ .* (Vp .- v)
    dθ = -v .* θ ./ dc .* @.log(v .* θ ./ dc)
    dv = 1 ./ (η ./ σn .+ a ./ v) .* (dτ ./ σn .- b .* dθ ./ θ)

    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dv # TODO: Flat fault only
    dvθ[2:3:end] .= 0 # TODO: Flat fault only
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_planarqdconst()
    # Constants and model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 500.00)
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    dc = 0.2
    blockvelx = 1e-9
    blockvely = 0.0

    # Create fault elements
    els = Elements()
    nfault = 1
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
    Vp = 1e-9
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * Vp / 1000 * ones(nnodes)
    # ics[2:3:end] = 1e-3 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)

    # Time integrate
    p = (∂t, els, η, dc)
    prob = ODEProblem(calc_dvθ1d, ics, tspan, p)
    sol = solve(prob, Rosenbrock23(autodiff = false), abstol = 1e-4, reltol = 1e-4)
    # sol = solve(prob, abstol = 1e-4, reltol = 1e-4)
    
    t = [x / siay for x in sol.t]
    θ = zeros(length(t), nfault)
    vx = zeros(length(t), nfault)
    vy = zeros(length(t), nfault)
    for i in 1:length(t)
        vx[i, :] = sol.u[i][1:3:end]
        vy[i, :] = sol.u[i][2:3:end]
        θ[i, :] = sol.u[i][3:3:end]
    end

    close("all")
    figure(figsize = (6, 12))
    subplot(2, 1, 1)
    plot(t, θ, "-b", linewidth = 0.5)
    xlabel("t (years)")
    ylabel(L"\theta")

    subplot(2, 1, 2)
    plot(t, vx, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("t (years)")
    ylabel(L"v_x")

    show()
end
ex_planarqdconst()
