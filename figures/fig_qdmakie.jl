using Revise
using DifferentialEquations
using AbstractPlotting
using Makie
using Random
using Printf
using NPZ
using Infiltrator
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
            plot!(subscene, [els.x1[idx[i]], els.x2[idx[i]]], [els.y1[idx[i]], els.y2[idx[i]]], color = :black)
        end
    end
    return nothing
end

function derivsconst(u, p, t)
    intidx, partials, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal = p
    nintidx = length(intidx)
    vxglobal = @. abs(u[1:3:end])
    vyglobal = u[2:3:end]
    theta = @. abs(u[3:3:end])
    dtracglobaldt =  partials["trac"]["fault"]["fault"] * [blockvxglobal .- vxglobal blockvyglobal .- vyglobal]'[:]
    vx, vy = multmatvec(els.rotmat[intidx, :, :], vxglobal, vyglobal)
    dtracxglobaldt, dtracyglobaldt = multmatvec(els.rotmat[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])

    dthetadt = zeros(nintidx)
    dvxdt = zeros(nintidx)
    dvydt = zeros(nintidx)
    for i in 1:length(intidx)
        dthetadt[i] = thetalaw(vx[i], theta[i], dc)
        dvxdt[i] = 1 / (eta / els.normalstress[intidx[i]] + els.a[intidx[i]] / vx[i]) * (dtracxglobaldt[i] / els.normalstress[intidx[i]] - els.b[intidx[i]] * dthetadt[i] / theta[i])
        dvydt[i] = 0
    end

    dvxglobaldt, dvyglobaldt = multmatvec(els.rotmatinv[intidx, :, :], dvxdt, dvydt)
    dudt = zeros(length(u))
    dudt[1:3:end] = dvxglobaldt
    dudt[2:3:end] = dvyglobaldt
    dudt[3:3:end] = dthetadt
    return dudt
end

function fig_qdmakie()
    # Fun things to play with
    nsteps = 500
    amplitude = 1.0
    nfault = 200

    # Constants
    siay = 365.25 * 24 * 60 * 60
    tspan = (0, siay * 100)
    abstol = 1e-4
    reltol = 1e-4
    mu = 3e10
    nu = 0.25
    eta = 0.25
    rho = 2700.0
    eta = mu / (2.0 * sqrt(mu / rho))
    dc = 0.05
    blockvelx = 1e-9
    blockvely = 0.0

    # Create fault elements
    els = Elements(Int(1e5))
    nnodes = 1 * nfault
    faultwidth = 10000
    x1, y1, x2, y2 = discretizedline(-faultwidth, 0, faultwidth, 0, nfault)

    # Spatially variable properties
    avec = collect(LinRange(0.010, 0.020, nfault))
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.a[els.endidx + i] = avec[i]
        els.b[els.endidx + i] = 0.020
        els.normalstress[els.endidx + i] = 50e6
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Convenience data structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    
    # Calculate slip to traction partials on the fault
    println("Calculating velocity to traction matrix")
    @time _, _, partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    
    # Set initial conditions
    ics = zeros(3 * nnodes)
    ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
    ics[2:3:end] = 0.0 * blockvely * ones(nnodes)
    ics[3:3:end] = 1e8 * ones(nnodes)
    intidx = collect(1:1:els.endidx) # indices of elements to integrate
    p = (intidx, partialsconst, els, eta, thetaaginglaw, dc, blockvelx, blockvely)
    prob = ODEProblem(derivsconst, ics, tspan, p)
    integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
    
    println("Step integrating")
    xplot = collect(1:1:nfault)
    textsize = 1.75
    vgloballimits = FRect(0, -2.0, nfault, 4.0)
    thetalimits = FRect(0, -1.0, nfault, 11.0)
    
    vxglobalupdate = Node(0.0 .* ones(nfault))
    vxglobalupdatemax = Node(0.0 .* ones(nfault))
    vxglobalupdatemaxvals = 0.0 .* ones(nfault)
    vyglobalupdate = Node(0.0 .* ones(nfault))
    vyglobalupdatemax = Node(0.0 .* ones(nfault))
    vyglobalupdatemaxvals = 0.0 .* ones(nfault)
    vxupdate = Node(0.0 .* ones(nfault))
    vxupdatemax = Node(0.0 .* ones(nfault))
    vxupdatemaxvals = 0.0 .* ones(nfault)
    vyupdate = Node(0.0 .* ones(nfault))
    vyupdatemax = Node(0.0 .* ones(nfault))
    vyupdatemaxvals = 0.0 .* ones(nfault)
    thetaupdate = Node(0.0 .* ones(nfault))
    thetaupdatemax = Node(0.0 .* ones(nfault))
    thetaupdatemaxvals = log10.(1e-10 .* ones(nfault))
    vglobalcurrent = Node(string(0))
    vcurrent = Node(string(0))
    thetacurrent = Node(string(0))

    scene = Scene(resolution=(3000, 2000))
    subscene1 = Scene(scene, IRect(50, 1050, 1350, 800))
    subscene2 = Scene(scene, IRect(1600, 1050, 1350, 800))
    subscene3 = Scene(scene, IRect(50, 50, 1350, 800))
    subscene4 = Scene(scene, IRect(1600, 50, 1350, 800))
    plotelementsmakie(els, subscene1)
    subscene1[Axis][:names][:axisnames] = ("x (m)", "y (m)")
    plot!(subscene2, xplot, vxglobalupdate, limits=vgloballimits, color = :red)
    plot!(subscene2, xplot, vxglobalupdatemax, limits=vgloballimits, color = :red, linestyle=:dot)
    plot!(subscene2, xplot, vyglobalupdate, limits=vgloballimits, color = :blue)
    plot!(subscene2, xplot, vyglobalupdatemax, limits=vgloballimits, color = :blue, linestyle=:dot)
    text!(subscene2, vglobalcurrent, position = (0.0, 2.0), align = (:left,  :center), textsize = textsize, limits=vgloballimits)
    subscene2[Axis][:names][:axisnames] = ("element index", "pow(v global)")
    plot!(subscene4, xplot, vxupdate, limits=vgloballimits, color = :red)
    plot!(subscene4, xplot, vxupdatemax, limits=vgloballimits, color = :red, linestyle=:dot)
    plot!(subscene4, xplot, vyupdate, limits=vgloballimits, color = :blue)
    plot!(subscene4, xplot, vyupdatemax, limits=vgloballimits, color = :blue, linestyle=:dot)
    text!(subscene4, vcurrent, position = (0.0, 2.0), align = (:left,  :center), textsize = textsize, limits=vgloballimits)
    subscene4[Axis][:names][:axisnames] = ("element index", "pow(v global)")
    plot!(subscene3, xplot, thetaupdate, limits=thetalimits, color = :green)
    plot!(subscene3, xplot, thetaupdatemax, limits=thetalimits, color = :green, linestyle=:dot)
    text!(subscene3, thetacurrent, position = (0.0, 10.0), align = (:left, :center), textsize = textsize, limits=thetalimits)
    subscene3[Axis][:names][:axisnames] = ("element index", "log10(θ)")
    display(scene)

    @time for i in 1:nsteps
        DifferentialEquations.step!(integrator)

        # Global velocities
        vxglobalupdate[] = vplot(integrator.u[1:3:end])
        vyglobalupdate[] = vplot(integrator.u[2:3:end])
        vxglobalupdatemaxvals = dropdims(findmax([integrator.u[1:3:end] vxglobalupdatemaxvals], dims=2)[1], dims=2)
        vxglobalupdatemax[] = vplot(vxglobalupdatemaxvals)
        vyglobalupdatemaxvals = dropdims(findmax([integrator.u[2:3:end] vyglobalupdatemaxvals], dims=2)[1], dims=2)
        vyglobalupdatemax[] = vplot(vyglobalupdatemaxvals)
        vglobalcurrent[] = @sprintf("t = %012.6f, n = %07d, max(vx) = %01.5f, min(vx) = %01.5f, max(vy) = %01.5f, min(vy) = %01.5f", integrator.t / siay, i, maximum(integrator.u[1:3:end]), minimum(integrator.u[1:3:end]), maximum(integrator.u[2:3:end]), minimum(integrator.u[2:3:end]))

        # Local velocities
        vx, vy = multmatvec(els.rotmat[1:els.endidx, :, :], integrator.u[1:3:end], integrator.u[2:3:end])
        vxupdate[] = vplot(vx)
        vyupdate[] = vplot(vy)
        vxupdatemaxvals = dropdims(findmax([integrator.u[1:3:end] vxupdatemaxvals], dims=2)[1], dims=2)
        vxupdatemax[] = vplot(vxupdatemaxvals)
        vyupdatemaxvals = dropdims(findmax([integrator.u[2:3:end] vyupdatemaxvals], dims=2)[1], dims=2)
        vyupdatemax[] = vplot(vyupdatemaxvals)
        vcurrent[] = @sprintf("t = %012.6f, n = %07d, max(vx) = %01.5f, min(vx) = %01.5f, max(vy) = %01.5f, min(vy) = %01.5f", integrator.t / siay, i, maximum(integrator.u[1:3:end]), minimum(integrator.u[1:3:end]), maximum(integrator.u[2:3:end]), minimum(integrator.u[2:3:end]))

        # State
        thetaupdate[] = log10.(integrator.u[3:3:end])
        thetaupdatemaxvals = dropdims(findmax([integrator.u[3:3:end] thetaupdatemaxvals], dims=2)[1], dims=2)
        thetaupdatemax[] = log10.(thetaupdatemaxvals)
        thetacurrent[] = @sprintf("t = %012.6f, n = %07d, max(θ) = %012.2f, min(θ) = %012.2f", integrator.t / siay, i, maximum(integrator.u[3:3:end]), minimum(integrator.u[3:3:end]))

        sleep(1e-10)
    end

end
fig_qdmakie()
