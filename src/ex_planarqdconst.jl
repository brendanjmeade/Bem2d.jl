using Revise
using DifferentialEquations
using PyCall
using PyPlot
using Bem2d

function plottimeseries(sol)
    siay = 365.25 * 24 * 60 * 60
    nfault = Int(size(sol)[1] / 3)
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
    figure(figsize = (15, 8))

    subplot(3, 2, 1)
    plot(t, vx, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_x")
    subplot(3, 2, 3)
    plot(t, vy, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_y")
    subplot(3, 2, 5)
    plot(t, θ, "-", linewidth = 0.5)
    yscale("log")
    xlabel("t (years)")
    ylabel(L"\theta")

    subplot(3, 2, 2)
    plot(1:1:length(t), vx, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_x")
    subplot(3, 2, 4)
    plot(1:1:length(t), vy, "-", linewidth = 0.5)
    yscale("log")
    ylabel(L"v_y")
    subplot(3, 2, 6)
    plot(1:1:length(t), θ, "-", linewidth = 0.5)
    yscale("log")
    xlabel("time step #")
    ylabel(L"\theta")

    figure(figsize = (15, 5))
    plotme = @.log10(vx')
    contourf(plotme, 20)
    colorbar()
    contour(plotme, 20, linewidths = 0.5, linestyles = "solid", colors = "k")
    xlabel("time step")
    ylabel("element index")
    show()
end

# Derivatives to feed to ODE integrator
function calc_dvθ(vθ, p, t)
    ∂t, els, η, dc, blockvx, blockvy = p
    vx = vθ[1:3:end]
    vy = vθ[2:3:end]
    θ = vθ[3:3:end]

    # Fault parallel and perpendicular velocity components
    blockvs = zeros(size(vx))
    blockvt = zeros(size(vx))
    for i in 1:els.endidx
        blockvs[i], blockvt[i] = els.rotmat[i, :, :] * [blockvx ; blockvy]
    end
    vs, vt = multmatvec(els.rotmat[1:els.endidx, :, :], vx, vy)    

    # Global velocites on fault from fault parallel and perpendicular components
    vxs, vys = multmatvec(els.rotmatinv[1:els.endidx, :, :], vs, zeros(size(vs)))
    vxt, vyt = multmatvec(els.rotmatinv[1:els.endidx, :, :], zeros(size(vt)), vt)
    
    # Change in tractions due to strike slip motion only
    dt = ∂t * [blockvx .- vxs blockvy .- vys]'[:] # interleaving

    # Component of tractions parallel to fault trace
    ts, tt = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])    
    
    # Fault perpendicular component???
    # Locked, vt = 0
    # Creep, vt = projected component of block rate (vxt, vyt)?
    # Viscous, vt ∼ tt (velocity proportional to traction)
    

    # Frictional slip for fault parallel traction
    dτ = ts # Shear stress on fault
    dθ = zeros(els.endidx)
    dv = zeros(els.endidx)
    for i in 1:els.endidx
        dθ[i] = -vs[i] * θ[i] / dc * log(vs[i] * θ[i] / dc) # slip law
        # dθ[i] = 1 - θ[i] * vmag[i] / dc # Aging law
        dv[i] = 1 / (η / els.σn[i] + els.a[i] / vs[i]) * (dτ[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
    end
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dv # TODO: Flat fault only
    dvθ[2:3:end] .= 0 # TODO: Flat fault only
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_planarqdconst()
    # TODO: Quadratic farfield partials
    # TODO: Quadratic coincident partials
    # TODO: Take only fault parallel velocity and traction
    # TODO: Option: Creeping tensile
    # TODO: Option: Tensile-slip rate scales with strike-slip rate
    # TODO: Option: Tensile slip is viscous, v~σ

    # Constants and model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0, siay * 500)
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
    els = Elements()
    nfault = 10
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = discretizedline(-faultwidth, 0, faultwidth, 0, nfault)
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
    println("Calculating velocity to traction matrix")
    srcidx = findall(x->x == "fault", els.name)
    obsidx = srcidx
    @time ∂u, ∂σ, ∂t = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)

    # Set initial conditions
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 1e-3 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)

    # Time integrate elastic model
    prob = ODEProblem(calc_dvθ, ics, tspan, (∂t, els, η, dc, blockvelx, blockvely))
    println("Time integrating")
    @time sol = solve(prob, RK4(), abstol = abstol, reltol = reltol)
    plottimeseries(sol)
end
ex_planarqdconst()
