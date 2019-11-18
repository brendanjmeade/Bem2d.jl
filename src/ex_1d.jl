using DifferentialEquations
using PyCall
using PyPlot
using Infiltrator

function calca(a, v)
    alow = 0.045
    ahighdiff = 0.0319 - 0.0450
    switchrange = 1e-2 # m/s
    a = alow + 0.5 * ahighdiff * (1 + tanh((v - 5e-3) / switchrange))
    return a
end

function calcdvθexp!(dvθ, vθ, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ = p
    θ = vθ[1]
    v = vθ[2]
    a = calca(a, v)
    dvθ[1] = -v * θ / dc * log(v * θ / dc)
    dvθ[2] = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * dvθ[1] / θ)
    return nothing
end

function calcdvθref!(dvθ, vθ, p, t)
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
    abstol = 1e-4
    reltol = 1e-4

    @infiltrate
    
    # Initial conditions
    ics = [1e8; Vp / 1000]
    
    # Time integrate
    p = (dc, η, σn, a, b, μ, Vp, L, ρ)
    prob1 = ODEProblem(calcdvθref!, ics, tspan, p)
    @time sol1 = solve(prob1, RK4(), abstol = abstol, reltol = reltol)
    prob2 = ODEProblem(calcdvθexp!, ics, tspan, p)
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
    yscale("log")
    ylabel(L"\theta")

    subplot(2, 2, 2)
    plot(1:1:length(t1), θ1, "-b", linewidth = 0.5)
    plot(1:1:length(t2), θ2, "-r", linewidth = 0.5)
    yscale("log")
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
