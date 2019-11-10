using DifferentialEquations
using BenchmarkTools
using PyCall
using PyPlot

function calcdvθ(vθ, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ = p
    θ = vθ[1]
    v = vθ[2]
    dθ = -v * θ / dc * log(v * θ / dc)
    dv = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * dθ / θ)
    dvθ = [dθ, dv]
    return dvθ
end

function calcdvθ!(dvθ, vθ, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ = p
    θ = vθ[1]
    v = vθ[2]
    dvθ[1] = -v * θ / dc * log(v * θ / dc)
    dvθ[2] = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * dvθ[1] / θ)
    return nothing
end

function ex_1d()
    # Model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 20000.0)
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
    abstol = 1e-10
    reltol = 1e-10
    
    # Initial conditions
    ics = [1e8; Vp / 1000]
    
    # Time integrate
    p = (dc, η, σn, a, b, μ, Vp, L, ρ)
    prob1 = ODEProblem(calcdvθ, ics, tspan, p)
    @time sol1 = solve(prob1, RK4(), abstol = abstol, reltol = reltol)
    prob2 = ODEProblem(calcdvθ!, ics, tspan, p)
    @time sol2 = solve(prob2, RK4(), abstol = abstol, reltol = reltol)
    
    t1 = [x / siay for x in sol1.t]
    θ1 = [x[1] for x in sol1.u]
    v1 = [x[2] for x in sol1.u]
    t2 = [x / siay for x in sol2.t]
    θ2 = [x[1] for x in sol2.u]
    v2 = [x[2] for x in sol2.u]

    close("all")
    figure(figsize = (12, 6))
    subplot(2, 2, 1)
    plot(t1, θ1, "-b", linewidth = 0.5)
    plot(t2, θ2, "-r", linewidth = 0.5)
    ylabel(L"\theta")

    subplot(2, 2, 2)
    plot(1:1:length(t1), θ1, "-b", linewidth = 0.5)
    plot(1:1:length(t2), θ2, "-r", linewidth = 0.5)
    ylabel(L"\theta")

    subplot(2, 2, 3)
    plot(t1, v1, "-b", linewidth = 0.5)
    plot(t2, v2, "-r", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("v (m/s)")

    subplot(2, 2, 4)
    plot(1:1:length(t1), v1, "-b", linewidth = 0.5)
    plot(1:1:length(t2), v2, "-r", linewidth = 0.5)
    yscale("log")
    xlabel("time (step)")
    ylabel("v (m/s)")

    return nothing
end
ex_1d()
