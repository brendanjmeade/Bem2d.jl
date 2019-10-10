using Revise
using DifferentialEquations
using PyCall
using PyPlot
using AbstractPlotting
using Makie
using Colors
using ColorSchemes
using JLD2
using Dates
using UUIDs
using Bem2d

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
        vpara = abs.(vpara) # Chris R. suggestion
        θ = abs.(θ)
        dθ[i] = -vpara[i] * θ[i] / dc * log(vpara[i] * θ[i] / dc) # slip law
        # dθ[i] = 1 - θ[i] * vpara[i] / dc # Aging law
        dvpara[i] = 1 / (η / els.σn[i] + els.a[i] / vpara[i]) * (dtpara[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
        dvperp[i] = 0 # fault perpendicular creep
        # dvperp[i] = dvpara[i] # fault perpendicular velocity proportional to fault parallel velocity
    end
    dvx, dvy = Bem2d.multmatvec(els.rotmatinv[1:els.endidx, :, :], dvpara, dvperp)
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dvx
    dvθ[2:3:end] = dvy
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_qdmakie()
    # Constants and model parameters
    PyPlot.close("all")
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
    els = Bem2d.Elements(Int(1e5))
    nfault = 100
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = Bem2d.discretizedline(-faultwidth, 0, faultwidth, 0, nfault)

    # Modify y1, and y2 for a sinusoidal fault
    amplitude = 1000.0
    y1 = amplitude * @.sin(2 * π * x1 / faultwidth)
    y2 = amplitude * @.sin(2 * π * x2 / faultwidth)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.a[els.endidx + 1:els.endidx + nfault] .= 0.015
    els.b[els.endidx + 1:els.endidx + nfault] .= 0.020
    els.σn[els.endidx + 1:els.endidx + nfault] .= 50e6
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    Bem2d.standardize_elements!(els)

    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    srcidx = findall(x->x == "fault", els.name)
    obsidx = srcidx
    @time ∂u, ∂σ, ∂t = Bem2d.∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)

    # Set initial conditions
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)

    p = (∂t, els, η, dc, blockvelx, blockvely)
    prob = DifferentialEquations.ODEProblem(calc_dvθ, ics, tspan, p)

    # (Step) Time integrate elastic model
    println("Step integrating")
    limits = Makie.FRect(0, -14, nfault, 20)
    xplot = collect(1:1:nfault)
    vxupdate = Makie.Node(1e-14 .* ones(nfault))
    vyupdate = Makie.Node(1e-14 .* ones(nfault))
    currenttimestep = Makie.Node(string(0))
    currenttime = Makie.Node(string(0))
    currentminvx = Makie.Node(string(0))
    currentmaxvx = Makie.Node(string(0))
    currentminvy = Makie.Node(string(0))
    currentmaxvy = Makie.Node(string(0))
    currentminθ = Makie.Node(string(0))
    currentmaxθ = Makie.Node(string(0))
    scene = Makie.plot(xplot, vxupdate, limits=limits, color = :red)
    Makie.plot!(xplot, vyupdate, limits=limits, color = :blue)
    Makie.text!(currenttime, position = (0, 6), align = (:left,  :center), textsize = 2, limits=limits)
    Makie.text!(currenttimestep, position = (0, 5), align = (:left,  :center), textsize = 2, limits=limits)
    Makie.text!(currentminvx, position = (nfault, 6), align = (:right,  :center), textsize = 2, limits=limits)
    Makie.text!(currentmaxvx, position = (nfault, 5), align = (:right,  :center), textsize = 2, limits=limits)
    Makie.text!(currentminvy, position = (nfault, 4), align = (:right,  :center), textsize = 2, limits=limits)
    Makie.text!(currentmaxvy, position = (nfault, 3), align = (:right,  :center), textsize = 2, limits=limits)
    Makie.text!(currentminθ, position = (nfault, 2), align = (:right,  :center), textsize = 2, limits=limits)
    Makie.text!(currentmaxθ, position = (nfault, 1), align = (:right,  :center), textsize = 2, limits=limits)
    axis = scene[Axis] # get the axis object from the scene
    axis[:names][:axisnames] = ("element index", "log v")
    Makie.display(scene)

    nsteps = 10000
    t = zeros(nsteps)
    vx = zeros(nsteps, nfault)
    vy = zeros(nsteps, nfault)
    θ = zeros(nsteps, nfault)
    integrator = init(prob, DP5(), abstol = abstol, reltol = reltol, progress = true)

    @time for i in 1:nsteps
        DifferentialEquations.step!(integrator)

        # Update Makie plot
        vxupdate[] = log10.(abs.(integrator.u[1:3:end]))
        vyupdate[] = log10.(abs.(integrator.u[2:3:end]))
        currenttime[] = string(integrator.t / siay)
        currenttimestep[] = string(i)
        currentminvx[] = string(minimum(integrator.u[1:3:end]))
        currentmaxvx[] = string(maximum(integrator.u[1:3:end]))
        currentminvy[] = string(minimum(integrator.u[2:3:end]))
        currentmaxvy[] = string(maximum(integrator.u[2:3:end]))
        currentminθ[] = string(minimum(integrator.u[3:3:end]))
        currentmaxθ[] = string(maximum(integrator.u[3:3:end]))
        sleep(1e-10)
    end
end
ex_qdmakie()
