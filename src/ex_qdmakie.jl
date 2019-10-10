using Revise
using DifferentialEquations
using AbstractPlotting
using Makie
using Printf
using Bem2d

function calc_dvθ(vθ, p, t)
    ∂t, els, η, dc, blockvx, blockvy = p
    vx = vθ[1:3:end]
    vy = vθ[2:3:end]
    θ = vθ[3:3:end]
    vpara, vperp = multmatvec(els.rotmat[1:els.endidx, :, :], vx, vy)
    dt = ∂t * [blockvx .- vx blockvy .- vy]'[:]
    dtpara, dtperp = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])
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
    # Fun things to play with
    nsteps = 1000
    amplitude = 0.0
    nfault = 100

    # Constants
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
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = Bem2d.discretizedline(-faultwidth, 0, faultwidth, 0, nfault)

    # Modify y1, and y2 for a sinusoidal fault
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
    integrator = init(prob, DP5(), abstol = abstol, reltol = reltol, progress = true)

    println("Step integrating")
    xplot = collect(1:1:nfault)
    textsize = 2
    vispower = 1.0/10.0
    limits = Makie.FRect(0, -2.0, nfault, 4.0)
    vxupdate = Makie.Node(0.0 .* ones(nfault))
    vxupdatemax = Makie.Node(0.0 .* ones(nfault))
    vxupdatemaxvals = 0.0 .* ones(nfault)
    vyupdate = Makie.Node(0 .* ones(nfault))
    vyupdatemax = Makie.Node(0.0 .* ones(nfault))
    vyupdatemaxvals = 0.0 .* ones(nfault)


    currenttime = Makie.Node(string(0))
    currentvx = Makie.Node(string(0))
    currentvy = Makie.Node(string(0))
    currentθ = Makie.Node(string(0))
    scene = Scene(resolution=(2000, 1000))
    Makie.plot!(xplot, vxupdate, limits=limits, color = :red)
    Makie.plot!(xplot, vxupdatemax, limits=limits, color = :red, linestyle=:dot)
    Makie.plot!(xplot, vyupdate, limits=limits, color = :blue)
    Makie.plot!(xplot, vyupdatemax, limits=limits, color = :blue, linestyle=:dot)

    Makie.text!(currenttime, position = (0.0, 2.0), align = (:left,  :center), textsize = textsize, limits=limits)
    Makie.text!(currentvx, position = (0.0, 1.8), align = (:left,  :center), textsize = textsize, limits=limits)
    Makie.text!(currentvy, position = (0.0, 1.6), align = (:left,  :center), textsize = 2, limits=limits)
    axis = scene[Axis] # get the axis object from the scene
    axis[:names][:axisnames] = ("element index", "v^p")
    Makie.display(scene)

    @time for i in 1:nsteps
        DifferentialEquations.step!(integrator)

        # Current velocities
        vxupdate[] = sign.(integrator.u[1:3:end]) .* (abs.(integrator.u[1:3:end])).^(vispower)
        vyupdate[] = sign.(integrator.u[2:3:end]) .* (abs.(integrator.u[2:3:end])).^(vispower)

        # Update maximum velocities
        vxupdatemaxvals = findmax([integrator.u[1:3:end] vxupdatemaxvals], dims=2)[1]
        vxupdatemaxvals = dropdims(vxupdatemaxvals, dims=2)
        vxupdatemax[] = sign.(vxupdatemaxvals) .* (abs.(vxupdatemaxvals)).^(vispower)
        vyupdatemaxvals = findmax([integrator.u[2:3:end] vyupdatemaxvals], dims=2)[1]
        vyupdatemaxvals = dropdims(vyupdatemaxvals, dims=2)
        vyupdatemax[] = sign.(vyupdatemaxvals) .* (abs.(vyupdatemaxvals)).^(vispower)

        # Update text
        currenttime[] = Printf.@sprintf("t = %012.6f, n = %d", integrator.t / siay, i)
        currentvx[] = Printf.@sprintf("max(vx) = %01.5f, min(vx) = %01.5f",
            maximum(integrator.u[1:3:end]), minimum(integrator.u[1:3:end]))
        currentvy[] = Printf.@sprintf("max(vy) = %01.5f, min(vy) = %01.5f",
            maximum(integrator.u[2:3:end]), minimum(integrator.u[2:3:end]))

        # currentminvx[] = string(minimum(integrator.u[1:3:end]))
        # currentmaxvx[] = string(maximum(integrator.u[1:3:end]))
        # currentminvy[] = string(minimum(integrator.u[2:3:end]))
        # currentmaxvy[] = string(maximum(integrator.u[2:3:end]))
        # currentminθ[] = string(minimum(integrator.u[3:3:end]))
        # currentmaxθ[] = string(maximum(integrator.u[3:3:end]))
        sleep(1e-10)
    end
end
ex_qdmakie()
