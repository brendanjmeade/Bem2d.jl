using Revise
using DifferentialEquations
using PyCall
using PyPlot

# Derivatives to feed to ODE integrator (returned values)
function calcdvθ(vθ, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ = p
    θ = vθ[1]
    v = vθ[2]
    dθ = -v * θ / dc * log(v * θ / dc)
    dv = 1 / (η / σn + a / v) * (μ * (Vp - v)  / (L * σn) - b * dθ / θ)
    dvθ = [dθ, dv]
    return dvθ
end

# Doesn't wor...yet: Derivatives to feed to ODE integrator (in place)
function calcdvθ!(dvθ, vθ, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ = p
    θ = vθ[1]
    v = vθ[2]
    dθ = -v * θ / dc * log(v * θ / dc)
    dv = 1 / (η / σn + a / v) * (μ * (Vp - v)  / (L * σn) - b * dθ / θ)
    dvθ = [dθ, dv]
    return nothing
end

function ex_1d()
    # Constants and model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 5000.0)
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    L = 60 * 1e3
    a = 0.015
    b = 0.02
    Vp = 1e-9
    σn = 30e6
    dc = 0.2

    # Set initial conditions
    ics = [1e8, Vp / 1000]

    # Time integrate
    p = (dc, η, σn, a, b, μ, Vp, L, ρ)
    prob = ODEProblem(calcdvθ, ics, tspan, p)
    sol = solve(prob, Tsit5(), abstol = 1e-4, reltol = 1e-4)    
    t = [x / siay for x in sol.t]
    θ = [x[1] for x in sol.u]
    v = [x[2] for x in sol.u]
    
    close("all")
    figure(figsize = (12, 6))
    subplot(2, 2, 1)
    plot(t, θ, "-b", linewidth = 0.5)
    ylabel(L"\theta")

    subplot(2, 2, 2)
    plot(1:1:length(t), θ, "-b", linewidth = 0.5)
    ylabel(L"\theta")

    subplot(2, 2, 3)
    plot(t, v, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("v (m/s)")

    subplot(2, 2, 4)
    plot(1:1:length(t), v, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (step)")
    ylabel("v (m/s)")
    return nothing
end
ex_1d()