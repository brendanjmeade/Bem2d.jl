using Revise
using DifferentialEquations
using AbstractPlotting
using Makie
using MakieThemes
using Printf
using Bem2d

function vplot(v)
    vispower = 1.0 / 10.0
    return sign.(v) .* (abs.(v)).^(vispower)
end

function plotelementsmakie(els, subscene)
    names = unique(els.name)
    for j in 1:length(names) - 1
        idx = findall(x->x == names[j], els.name)
        for i in 1:length(idx)
            Makie.plot!(subscene, [els.x1[idx[i]], els.x2[idx[i]]], [els.y1[idx[i]], els.y2[idx[i]]], color = :black)
        end
    end
    return nothing
end

function calc_dvθ(vθ, p, t)
    ∂t, els, η, dc, blockvxglobal, blockvyglobal = p
    vxglobal = vθ[1:3:end]
    vyglobal = vθ[2:3:end]
    θ = vθ[3:3:end]
    vx, vy = multmatvec(els.rotmat[1:els.endidx, :, :], vxglobal, vyglobal)
    dt = ∂t * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]
    dtx, dty = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])
    dθ, dvx, dvy = zeros(els.endidx), zeros(els.endidx), zeros(els.endidx)
    for i in 1:els.endidx
        vx = abs.(vx) # Chris R. suggestion
        θ = abs.(θ) # Chris R. suggestion
        dθ[i] = -vx[i] * θ[i] / dc * log(vx[i] * θ[i] / dc) # slip law
        # dθ[i] = 1 - θ[i] * vx[i] / dc # Aging law
        dvx[i] = 1 / (η / els.σn[i] + els.a[i] / vx[i]) * (dtx[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
        dvy[i] = 0 # fault perpendicular creep
        # dvperp[i] = dvpara[i] # fault perpendicular velocity proportional to fault parallel velocity
    end
    dvxglobal, dvyglobal = Bem2d.multmatvec(els.rotmatinv[1:els.endidx, :, :], dvx, dvy)
    dvθ = zeros(3 * els.endidx)
    dvθ[1:3:end] = dvxglobal
    dvθ[2:3:end] = dvyglobal
    dvθ[3:3:end] = dθ
    return dvθ
end

function ex_qdmakie()
    # Fun things to play with
    nsteps = 5000
    amplitude = 100.0
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
    y1 = amplitude * @.sin(8 * π * x1 / faultwidth)
    y2 = amplitude * @.sin(8 * π * x2 / faultwidth)
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
    # integrator = init(prob, DP5(), abstol = abstol, reltol = reltol, progress = true)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol, progress = true)

    println("Step integrating")
    xplot = collect(1:1:nfault)
    textsize = 1.75
    vgloballimits = Makie.FRect(0, -2.0, nfault, 4.0)
    θlimits = Makie.FRect(0, -1.0, nfault, 11.0)

    vxglobalupdate = Makie.Node(0.0 .* ones(nfault))
    vxglobalupdatemax = Makie.Node(0.0 .* ones(nfault))
    vxglobalupdatemaxvals = 0.0 .* ones(nfault)
    vyglobalupdate = Makie.Node(0.0 .* ones(nfault))
    vyglobalupdatemax = Makie.Node(0.0 .* ones(nfault))
    vyglobalupdatemaxvals = 0.0 .* ones(nfault)

    vxupdate = Makie.Node(0.0 .* ones(nfault))
    vxupdatemax = Makie.Node(0.0 .* ones(nfault))
    vxupdatemaxvals = 0.0 .* ones(nfault)
    vyupdate = Makie.Node(0.0 .* ones(nfault))
    vyupdatemax = Makie.Node(0.0 .* ones(nfault))
    vyupdatemaxvals = 0.0 .* ones(nfault)

    θupdate = Makie.Node(0.0 .* ones(nfault))
    θupdatemax = Makie.Node(0.0 .* ones(nfault))
    θupdatemaxvals = log10.(1e-10 .* ones(nfault))

    vglobalcurrent = Makie.Node(string(0))
    vcurrent = Makie.Node(string(0))
    θcurrent = Makie.Node(string(0))

    scene = Makie.Scene(resolution=(3000, 2000))
    subscene1 = Makie.Scene(scene, Makie.IRect(50, 1050, 1350, 800))
    subscene2 = Makie.Scene(scene, Makie.IRect(1600, 1050, 1350, 800))
    subscene3 = Makie.Scene(scene, Makie.IRect(50, 50, 1350, 800))
    subscene4 = Makie.Scene(scene, Makie.IRect(1600, 50, 1350, 800))

    plotelementsmakie(els, subscene1)
    subscene1[Axis][:names][:axisnames] = ("x (m)", "y (m)")

    Makie.plot!(subscene2, xplot, vxglobalupdate, limits=vgloballimits, color = :red)
    Makie.plot!(subscene2, xplot, vxglobalupdatemax, limits=vgloballimits, color = :red, linestyle=:dot)
    Makie.plot!(subscene2, xplot, vyglobalupdate, limits=vgloballimits, color = :blue)
    Makie.plot!(subscene2, xplot, vyglobalupdatemax, limits=vgloballimits, color = :blue, linestyle=:dot)
    Makie.text!(subscene2, vglobalcurrent, position = (0.0, 2.0), align = (:left,  :center), textsize = textsize, limits=vgloballimits)
    subscene2[Axis][:names][:axisnames] = ("element index", "pow(v global)")

    Makie.plot!(subscene4, xplot, vxupdate, limits=vgloballimits, color = :red)
    Makie.plot!(subscene4, xplot, vxupdatemax, limits=vgloballimits, color = :red, linestyle=:dot)
    Makie.plot!(subscene4, xplot, vyupdate, limits=vgloballimits, color = :blue)
    Makie.plot!(subscene4, xplot, vyupdatemax, limits=vgloballimits, color = :blue, linestyle=:dot)
    Makie.text!(subscene4, vcurrent, position = (0.0, 2.0), align = (:left,  :center), textsize = textsize, limits=vgloballimits)
    subscene4[Axis][:names][:axisnames] = ("element index", "pow(v global)")

    Makie.plot!(subscene3, xplot, θupdate, limits=θlimits, color = :green)
    Makie.plot!(subscene3, xplot, θupdatemax, limits=θlimits, color = :green, linestyle=:dot)
    Makie.text!(subscene3, θcurrent, position = (0.0, 10.0), align = (:left,  :center), textsize = textsize, limits=θlimits)
    subscene3[Axis][:names][:axisnames] = ("element index", "log10(θ)")
    Makie.display(scene)


    @time for i in 1:nsteps
        DifferentialEquations.step!(integrator)

        # Global velocities
        vxglobalupdate[] = vplot(integrator.u[1:3:end])
        vyglobalupdate[] = vplot(integrator.u[2:3:end])
        vxglobalupdatemaxvals = dropdims(findmax([integrator.u[1:3:end] vxglobalupdatemaxvals], dims=2)[1], dims=2)
        vxglobalupdatemax[] = vplot(vxglobalupdatemaxvals)
        vyglobalupdatemaxvals = dropdims(findmax([integrator.u[2:3:end] vyglobalupdatemaxvals], dims=2)[1], dims=2)
        vyglobalupdatemax[] = vplot(vyglobalupdatemaxvals)
        vglobalcurrent[] = Printf.@sprintf("t = %012.6f, n = %07d, max(vx) = %01.5f, min(vx) = %01.5f, max(vy) = %01.5f, min(vy) = %01.5f", integrator.t / siay, i, maximum(integrator.u[1:3:end]), minimum(integrator.u[1:3:end]), maximum(integrator.u[2:3:end]), minimum(integrator.u[2:3:end]))

        # Local velocities
        vx, vy = Bem2d.multmatvec(els.rotmat[1:els.endidx, :, :], integrator.u[1:3:end], integrator.u[2:3:end])
        vxupdate[] = vplot(vx)
        vyupdate[] = vplot(vy)
        vxupdatemaxvals = dropdims(findmax([integrator.u[1:3:end] vxupdatemaxvals], dims=2)[1], dims=2)
        vxupdatemax[] = vplot(vxupdatemaxvals)
        vyupdatemaxvals = dropdims(findmax([integrator.u[2:3:end] vyupdatemaxvals], dims=2)[1], dims=2)
        vyupdatemax[] = vplot(vyupdatemaxvals)
        vcurrent[] = Printf.@sprintf("t = %012.6f, n = %07d, max(vx) = %01.5f, min(vx) = %01.5f, max(vy) = %01.5f, min(vy) = %01.5f", integrator.t / siay, i, maximum(integrator.u[1:3:end]), minimum(integrator.u[1:3:end]), maximum(integrator.u[2:3:end]), minimum(integrator.u[2:3:end]))



        # State
        θupdate[] = log10.(integrator.u[3:3:end])
        θupdatemaxvals = dropdims(findmax([integrator.u[3:3:end] θupdatemaxvals], dims=2)[1], dims=2)
        θupdatemax[] = log10.(θupdatemaxvals)
        θcurrent[] = Printf.@sprintf("t = %012.6f, n = %07d, max(θ) = %012.2f, min(θ) = %012.2f", integrator.t / siay, i, maximum(integrator.u[3:3:end]), minimum(integrator.u[3:3:end]))

        sleep(1e-10)
    end
end
ex_qdmakie()
