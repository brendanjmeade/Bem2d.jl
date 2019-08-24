using Revise
using PyCall
using PyPlot
using Bem2d

# Derivatives to feed to ODE integrator
function derivs(t, vθ, ∂t, els, blockvelx, blockvely)
    # Current shear stress rate on fault (slip -> traction)
    dτ = ∂t * [vθ[1:3:end] ; vθ[2:3:end]]
    dθ = 1 - θ * v / Dc # state evolution: aging law
    dθ = 1 - θ * v / Dc # state evolution: slip law
    dv = 1 / (η / els.σn[1:els.endidx] + els.a[1:els.endidx] / v) *
        (dτ / els.σn[1:els.endidx] - els.b[1:els.endidx] * dθ / θ) # velocity evolution
    
    derivs = zeros(3 * els.endidx)
    derivs[1:3:end] .= -currentvels[1:2:end] + blockvelx
    derivs[2:3:end] .= -currentvels[2:2:end] + blockvely
    derivs[3:3:end] .= dθ
    return derivs
end

function ex_planarqdconst()
    # Constants and model parameters
    spy = 365.25 * 24 * 60 * 60  # Seconds per year
    timeinterval = spy * [0.0 1000.0]
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    Dc = 0.05
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
    return nothing
end
ex_planarqdconst()

# function plottimeseries(solution, spy)
#     figure(figsize=(12, 9))
#     subplot(3, 2, 1)
#     plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 0::3], linewidth=0.5)
#     ylabel(L"$u_x$ (m)")

#     subplot(3, 2, 2)
#     plot(SOLUTION["y"][:, 0::3], linewidth=0.5)
#     ylabel(L"$u_x$ (m)")

#     subplot(3, 2, 3)
#     plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 1::3], linewidth=0.5)
#     ylabel(L"$u_y$ (m)")

#     subplot(3, 2, 4)
#     plot(SOLUTION["y"][:, 1::3], linewidth=0.5)
#     ylabel(L"$u_y$ (m)")

#     subplot(3, 2, 5)
#     plot(SOLUTION["t"] / SPY, SOLUTION["y"][:, 2::3], linewidth=0.5)
#     xlabel("time (years)")
#     ylabel(L"$\theta$")

#     subplot(3, 2, 6)
#     plot(SOLUTION["y"][:, 2::3], linewidth=0.5)
#     xlabel("steps")
#     ylabel(L"$\theta$")
#     show()
# end
# plottimeseries()