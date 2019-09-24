using Revise
using DifferentialEquations
using AbstractPlotting
using Makie
using Bem2d

# Derivatives to feed to ODE integrator
function calc_dvθ(vθ, p, t)
    ∂t, els, η, dc, blockvx, blockvy = p
    vx = vθ[1:3:end]
    vy = vθ[2:3:end]
    θ = vθ[3:3:end]
    println("In derivative calcuation")
    @show vx

    # Global fault velocites from fault parallel and perpendicular components
    vpara, vperp = multmatvec(els.rotmat[1:els.endidx, :, :], vx, vy)

    # Change in tractions (global then local)
    dt = ∂t * [blockvx .- vx blockvy .- vy]'[:]
    dtpara, dtperp = multmatvec(els.rotmat[1:els.endidx, :, :], dt[1:2:end], dt[2:2:end])

    # Frictional slip for fault parallel traction
    dθ, dvpara, dvperp = zeros(els.endidx), zeros(els.endidx), zeros(els.endidx)

    for i in 1:els.endidx
        dθ[i] = -vpara[i] * θ[i] / dc * log(vpara[i] * θ[i] / dc) # slip law
        dvpara[i] = 1 / (η / els.σn[i] + els.a[i] / vpara[i]) * (dtpara[i] / els.σn[i] - els.b[i] * dθ[i] / θ[i])
        dvperp[i] = 0 # fault perpendicular creep
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
    abstol = 1e-6
    reltol = 1e-6
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
    prob = ODEProblem(calc_dvθ, ics, tspan, p)

    # (Step) Time integrate elastic model
    println("Step integrating")

    # Set up Makie plot for realtime visualization
    limits = FRect(0, -14, nfault, 20)
    xplot = collect(1:1:nfault)
    vxupdate = Node(1e-14 .* ones(nfault))
    currenttimestep = Node(string(0))
    currenttime = Node(string(0))
    currentminvx = Node(string(0))
    currentmaxvx = Node(string(0))
    currentminvy = Node(string(0))
    currentmaxvy = Node(string(0))
    currentminθ = Node(string(0))
    currentmaxθ = Node(string(0))
    scene = Makie.plot(xplot, vxupdate, limits=limits, color = :red)
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
    display(scene)

    ###
    ### The RK4 integrator works well
    ### DP8 and a bunch of others fail while iterating in the derivatives
    ### function.  They fail because they produce an intermediate velocity
    ### that is negative.  This is a problem because the derivatives require
    ### the log of this value.
    ###
    ### Help!
    ###
    integrator = init(prob, RK4(), abstol = abstol, reltol = reltol, progress = true)
    # integrator = init(prob, VCABM5(), abstol = abstol, reltol = reltol, progress = true)
    # integrator = init(prob, DP5(), abstol = abstol, reltol = reltol, progress = true)
    # integrator = init(prob, DP8(), abstol = abstol, reltol = reltol, progress = true)
    nsteps = 1000

    @time for i in 1:nsteps
        println("In for loop before integrator step")
        @show i
        @show integrator.u[1:3:end]
        DifferentialEquations.step!(integrator)

        # Update Makie plot
        vxupdate[] = log10.(abs.(integrator.u[1:3:end]))
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
ex_planarqdconst()
