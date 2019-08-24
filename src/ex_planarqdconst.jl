using Revise
using NLsolve
using PyCall
using PyPlot
using Roots
using Bem2d

# State evolution law: aging law
function calcdθdt!(dθdt, θ, v, b, v0, Dc, f0)
    dθdt = b .* v0 ./ Dc .* (@.exp((f0 .- θ) ./ b) .- (v ./ v0))
end

# https://discourse.julialang.org/t/parametrize-nlsolve-function/11886

# Should I eliminate this since its only called once?
# Steady-state state for initial condition """
# function steadystate(vels)
#     θ = ones(size(vels))
#     res = nlsolve((F, x) -> f!(F, x, a, b), [0.0], autodiff=:forward).zero
#     return θ
# end

# Find the sliding velocity that balances 
function f(v, τ, σn, θ, η, a, v0)
    return τ - η * v - (σn * a * asinh(v / (2 * v0) * exp(θ / a)))
end


# Derivatives to feed to ODE integrator
function derivatives(t, xθ, ∂t, oldvels, els, nnodes, blockvelx, blockvely)
    # Current shear stress on fault (slip -> traction)
    t = ∂t * [xθ[1:3:end] ; xθ[2:3:end]]

    # Solve for the current velocity...This is the algebraic part
    currentvels = zeros(2 * nnodes)
    for i in 1:nnodes
        vels = fsolve(f, oldvels[2 * i - 1], args = (t[2 * i - 1], els.σn[i], xθ[3:3:end][i]))
        currentvels[2 * i - 1] = vels
        currentvels[2 * i] = 0  # Assuming x-direction only
    end

    dθdt = calcdθdt!(xθ[3:3:end], currentvels[1:2:end])
    derivs = zeros(3 * nnodes)
    derivs[1:3:end] = -currentvels[1:2:end] + blockvelx
    derivs[2:3:end] = -currentvels[2:2:end] + blockvely
    derivs[3:3:end] = dθdt
    return derivs
end

function ex_planarqdconst()
    # Constants and model parameters
    spy = 365.25 * 24 * 60 * 60  # Seconds per year
    timeinterval = spy * [0.0, 1000.0]
    μ = 3e10
    ν = 0.25
    ρ = 2700.0  # rock density (kg/m^3)
    η = μ / (2.0 * sqrt(μ / ρ)) # Radiation damping coefficient (kg / (m^2 * s))
    Dc = 0.05  # state evolution length scale (m)
    f0 = 0.6  # baseline coefficient of friction
    blockvelx = 1e-9
    blockvely = 0.0
    v0 = 1e-6  # reference velocity

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
    els.σn[els.endidx + 1:els.endidx + nfault] .= 50e6
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Calculate slip to traction partials on the fault
    println("Calculating ∂s")
    srcidx = findall(x->x == "fault", els.name)
    obsidx = srcidx
    _, _, ∂t = ∂constslip(els, srcidx, obsidx, μ, ν)

    # Try calculating derivatives that we want to integrate

    # Set initial conditions and time integrate
    initconds = zeros(3 * nnodes)
    initconds[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    initconds[2:3:end] = 1e-3 * blockvely * ones(nnodes)
    initconds[3:3:end] = 0.5 * ones(nnodes)

    # println("Solving for steady state initial conditions...no idea what Im doing")
    # res = nlsolve((F, x) -> dθdt!(F, x, initvelmag, els.b[1:els.endidx], v0, Dc, f0), initconds[3:3:end], autodiff=:forward).zero
    # res = nlsolve((dθdt, θ)->calcdθdt!(dθdt, θ, initvelmag, els.b[1:els.endidx], v0, Dc, f0), initconds[3:3:end], autodiff = :forward)
    # display(res)

end
ex_planarqdconst()

function plottimeseries(solution, spy)
    figure(figsize=(12, 9))
    subplot(3, 2, 1)
    plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 0::3], linewidth=0.5)
    ylabel(L"$u_x$ (m)")

    subplot(3, 2, 2)
    plot(SOLUTION["y"][:, 0::3], linewidth=0.5)
    ylabel(L"$u_x$ (m)")

    subplot(3, 2, 3)
    plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 1::3], linewidth=0.5)
    ylabel(L"$u_y$ (m)")

    subplot(3, 2, 4)
    plot(SOLUTION["y"][:, 1::3], linewidth=0.5)
    ylabel(L"$u_y$ (m)")

    subplot(3, 2, 5)
    plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 2::3], linewidth=0.5)
    xlabel("time (years)")
    ylabel(L"$\theta$")

    subplot(3, 2, 6)
    plot(SOLUTION["y"][:, 2::3], linewidth=0.5)
    xlabel("steps")
    ylabel(L"$\theta$")
    show()
end
plottimeseries()