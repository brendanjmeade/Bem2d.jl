using Revise
using DifferentialEquations
using PyCall
using PyPlot
using AbstractPlotting
using Makie

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
    contourf(plotme, 20, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}v$ (m/s)")
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

    # # Global block velocities from fault parallel and perpendicular velocity components
    # blockvpara, blockvperp = multmatvecsingle(els.rotmat[1:els.endidx, :, :], blockvx, blockvy)
    # blockvxpara, blockvypara = multmatvec(els.rotmatinv[1:els.endidx, :, :], blockvpara, zeros(size(blockvperp)))
    # blockvxperp, blockvyperp = multmatvec(els.rotmatinv[1:els.endidx, :, :], zeros(size(blockvpara)), blockvperp)

    # Global fault velocites from fault parallel and perpendicular components
    vpara, vperp = multmatvec(els.rotmat[1:els.endidx, :, :], vx, vy)
    # vxpara, vypara = multmatvec(els.rotmatinv[1:els.endidx, :, :], vpara, zeros(size(vperp)))
    # vxperp, vyperp = multmatvec(els.rotmatinv[1:els.endidx, :, :], zeros(size(vpara)), vperp)

    # Change in tractions (global then local)
    dt = ∂t * [blockvx .- vx blockvy .- vy]'[:]
    dtpara, dtperp = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])

    # Frictional slip for fault parallel traction
    dθ, dvpara, dvperp = zeros(els.endidx), zeros(els.endidx), zeros(els.endidx)
    for i in 1:els.endidx
        dθ[i] = -vpara[i] * θ[i] / dc * log(vpara[i] * θ[i] / dc) # slip law
        # dθ[i] = 1 - θ[i] * vpara[i] / dc # Aging law
        dvpara[i] = 1 / (η / els.σn[i] + els.a[i] / vpara[i]) * (dtpara[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
        dvperp[i] = 0 # fault perpendicular creep
        # dvperp[i] = dvpara[i] # fault perpendicular velocity proportional to fault parallel velocity
    end
    dvx, dvy = multmatvec(els.rotmatinv[1:els.endidx, :, :], dvpara, dvperp)
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dvx
    dvθ[2:3:end] = dvy
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_planarqdconst()
    # Constants and model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0, siay * 1000)
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
    nfault = 100
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = discretizedline(-faultwidth, 0, faultwidth, 0, nfault)

    # Modify y1, and y2 for a sinusoidal fault
    amplitude = 1000.0
    # y1 = amplitude * @.sin(2 * π * x1 / faultwidth)
    # y2 = amplitude * @.sin(2 * π * x2 / faultwidth)
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
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)

    # (Bulk) Time integrate elastic model
    p = (∂t, els, η, dc, blockvelx, blockvely)
    calc_dvθ(ics, p, 1)
    prob = ODEProblem(calc_dvθ, ics, tspan, p)
    println("Bulk integrating")
    # @time sol = solve(prob, RK4(), abstol = abstol, reltol = reltol, progress = true)
    # plottimeseries(sol)
#
    # TODO: These don't seem to work for non-planar faults
    # TODO: Try continuity to see if that makes a difference
    # @time sol = solve(prob, VCABM5(), abstol = abstol, reltol = reltol, progress = true)
    # @time sol = solve(prob, DP5(), abstol = abstol, reltol = reltol, progress = true)
    # @time sol = solve(prob, DP8(), abstol = abstol, reltol = reltol, progress = true)

    # (Step) Time integrate elastic model
    println("Step integrating")
    # Set up Makie plot for realtime visualization
    limits = FRect(0, -14, nfault, 20)
    xplot = collect(1:1:nfault)
    yupdate = Node(1e-14 .* ones(nfault))
    scene = Makie.plot(xplot, yupdate, limits=limits, color = :red)
    display(scene)

    # nsteps = length(sol.t)
    nsteps = 10000
    t = zeros(nsteps)
    vx = zeros(nsteps, nfault)
    vy = zeros(nsteps, nfault)
    θ = zeros(nsteps, nfault)
    integrator = init(prob, RK4(), abstol = abstol, reltol = reltol, progress = true)

    @time for i in 1:nsteps
        step!(integrator)
        t[i] = integrator.t
        vx[i, :] = integrator.u[1:3:end]
        vy[i, :] = integrator.u[2:3:end]
        θ[i, :] = integrator.u[3:3:end]
        yupdate[] = log10.(abs.(vx[i, :]))
        println(i)
    end

end
ex_planarqdconst()
