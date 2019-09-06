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

    # Global block velocities from fault parallel and perpendicular velocity components
    blockvpara, blockvperp = multmatvecsingle(els.rotmat[1:els.endidx, :, :], blockvx, blockvy)
    blockvxpara, blockvypara = multmatvec(els.rotmatinv[1:els.endidx, :, :], blockvpara, zeros(size(blockvperp)))
    blockvxperp, blockvyperp = multmatvec(els.rotmatinv[1:els.endidx, :, :], zeros(size(blockvpara)), blockvperp)

    # Global fault velocites from fault parallel and perpendicular components
    vpara, vperp = multmatvec(els.rotmat[1:els.endidx, :, :], vx, vy)
    vxpara, vypara = multmatvec(els.rotmatinv[1:els.endidx, :, :], vpara, zeros(size(vperp)))
    vxperp, vyperp = multmatvec(els.rotmatinv[1:els.endidx, :, :], zeros(size(vpara)), vperp)
    
    # Change in tractions due to fault parallel motion
    Δvx = blockvxpara - vxpara
    Δvy = blockvypara - vypara
    dtpara = ∂t * [Δvx Δvy]'[:]
    tparapara, tparaperp = multmatvec(els.rotmat[1:els.endidx, :, :], dtpara[1:2:end], dtpara[2:2:end])

    # # NO TENSILE SLIP - Change in tractions due to fault perpendicular motion
    # tperppara, tperpperp = zeros(size(vxperp)), zero(size(vxperp))

    # CREEP - Change in tractions due to fault perpendicular motion
    Δvx = blockvxperp - vxperp
    Δvy = blockvyperp - vyperp
    dtperp = ∂t * [Δvx Δvy]'[:]
    tperppara, tperpperp = multmatvec(els.rotmat[1:els.endidx, :, :], dtperp[1:2:end], dtperp[2:2:end])

    # # VISCOUS - Change in tractions due to fault perpendicular motion
    # viscosity = 1e19 # Pa⋅s
    # fault_width = 1.0 # meter
    # Δvx = blockvxperp - vxperp
    # Δvy = blockvyperp - vyperp
    # dtperp = ∂t * [Δvx Δvy]'[:]
    # tperppara, tperpperp = multmatvec(els.rotmat[1:els.endidx, :, :], dtperp[1:2:end], dtperp[2:2:end])

    # TOTAL - Change in tractions
    tpara = tparapara + tperppara
    tperp = tparaperp + tperpperp
        
    # Frictional slip for fault parallel traction
    dτ = tpara # Shear stress on fault
    dθ = zeros(els.endidx)
    dv = zeros(els.endidx)
    for i in 1:els.endidx
        dθ[i] = -vpara[i] * θ[i] / dc * log(vpara[i] * θ[i] / dc) # slip law
        # dθ[i] = 1 - θ[i] * vmag[i] / dc # Aging law
        dv[i] = 1 / (η / els.σn[i] + els.a[i] / vpara[i]) * (dτ[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
    end
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dv # TODO: Flat fault only
    dvθ[2:3:end] .= 0 # TODO: Flat fault only
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_planarqdconst()
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
    els = Elements(Int(1e5))
    nfault = 30
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
    p = (∂t, els, η, dc, blockvelx, blockvely)
    calc_dvθ(ics, p, 1)
    prob = ODEProblem(calc_dvθ, ics, tspan, p)
    println("Time integrating")
    @time sol = solve(prob, RK4(), abstol = abstol, reltol = reltol, progress = true)
    plottimeseries(sol)
end
ex_planarqdconst()
