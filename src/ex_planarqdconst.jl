using Revise
using DifferentialEquations
using PyCall
using PyPlot
using Bem2d

function plottimeseries(sol, nfault, siay)
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

# Derivatives to feed to ODE integrator
function calc_dvθ(vθ, p, t)
    display(t / (365.25 * 24 * 60 * 60))
    ∂t, els, η, dc = p
    θ = vθ[3:3:end]
    v = @.sqrt(vθ[1:3:end].^2 + vθ[2:3:end].^2) # TODO: Flat fault only
    dt = ∂t * [vθ[1:3:end] ; vθ[2:3:end]]
    dτ = dt[1:2:end]

    display(dτ)

    dθ = zeros(els.endidx)
    dv = zeros(els.endidx)
    for i in 1:els.endidx
        dθ[i] = -v[i] * θ[i] / dc * log(v[i] * θ[i] / dc)
        dv[i] = 1 / (η / els.σn[i] + els.a[i] / v[i]) * (dτ[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
    end

    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dv # TODO: Flat fault only
    dvθ[2:3:end] .= 0 # TODO: Flat fault only
    dvθ[3:3:end] = dθ
    return dvθ
end


# # 1-D Derivatives to feed to ODE integrator
# function calc_dvθ1d(vθ, p, t)
#     # display(t / (365.25 * 24 * 60 * 60))
#     v = vθ[1:3:end]
#     θ = vθ[3:3:end]
#     dc, η, μ, blockvelx, blockvely, L, ρ, els = p
#     siay = 365.25 * 24 * 60 * 60

#     # Hacky "1d" version
#     dθ = zeros(els.endidx)
#     dv = zeros(els.endidx)
#     for i in 1:els.endidx
#         dτ = μ * (blockvelx - v[i]) / L
#         dθ[i] = -v[i] * θ[i] / dc * log(v[i] * θ[i] / dc)
#         dv[i] = 1 / (η / els.σn[i] + els.a[i] / v[i]) * (dτ / els.σn[i] - els.b[i] * dθ[i] / θ[i])
#     end
    
#     dvθ = zeros(3 * els.endidx)
#     dvθ[1:3:end] = dv # TODO: Flat fault only
#     dvθ[2:3:end] .= 0 # TODO: Flat fault only
#     dvθ[3:3:end] = dθ
#     return dvθ
# end

function ex_planarqdconst()
    # Constants and model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 505.00)
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    dc = 0.2
    blockvelx = 1e-9
    blockvely = 0.0
    L = 60 * 1e3
    Vp = 1e-9

    # Create fault elements
    els = Elements()
    nfault = 3
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = discretizedline(-faultwidth, 0, faultwidth, 0, nfault)
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
    ∂u, ∂σ, ∂t = ∂constslip(els, srcidx, obsidx, μ, ν)

    display(∂σ)  # Rows look strange here
    close("all")
    matshow(∂u)
    colorbar()

    matshow(∂σ)
    colorbar()

    matshow(∂t)
    colorbar()
    show()
    # Set initial conditions
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 1e-3 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)

    # # Time integrate psuedo 1d model
    # p1d = (dc, η, μ, blockvelx, blockvely, L, ρ, els)
    # prob = ODEProblem(calc_dvθ1d, ics, tspan, p1d)
    # @time sol = solve(prob, RK4(), abstol = 1e-4, reltol = 1e-4)



    # Time integrate elastic model
    # p = (∂t, els, η, dc)
    # prob = ODEProblem(calc_dvθ, ics, tspan, p)
    # @time sol = solve(prob, RK4(), abstol = 1e-4, reltol = 1e-4)
    
    # plottimeseries(sol, nfault, siay)
end
ex_planarqdconst()
