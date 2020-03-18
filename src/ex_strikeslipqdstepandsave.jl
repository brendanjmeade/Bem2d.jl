using Revise
using DifferentialEquations
using JLD2
using Dates
using Infiltrator
using PyPlot

function stressfunction(slip, mu, x, y, d1, d2)
    stressd1 = @. slip*mu/(2*pi) * ((y+d1)/(x^2+(y+d1)^2) - (y+d1)/(x^2+(y-d1)^2))
    stressd2 = @. slip*mu/(2*pi) * ((y+d2)/(x^2+(y+d2)^2) - (y+d2)/(x^2+(y-d2)^2))
    stress = stressd1 .- stressd2
    return stress
end

function interactionmatrix(nels, dtop, dbot, mu)
    mat = zeros(nels, nels) # to store results
    dvec = collect(LinRange(dtop, dbot, nels+1))
    dtopvec = dvec[1:1:end-1]
    dbotvec = dvec[2:1:end]
    ymid = @. -(dtopvec + dbotvec) / 2 
    
    # Loop over all possible combinations
    for i = 1:nels
        for j = 1:nels
            mat[i, j] = stressfunction(1, mu, 0, ymid[i], dtopvec[j], dbotvec[j])
        end
    end
    return mat
end

function strikeslipstress()
    PyPlot.close("all")
    nels = 100
    mu = 3e10
    d1 = 20e3 # meters
    d2 = 10e3 # meters

    #! Calculate full element to element interaction interaction
    @time mat = interactionmatrix(nels, 0e3, 20e3, mu)
    PyPlot.matshow(mat) # Plot interaction matrix
end
strikeslipstress()


function derivsconst!(dudt, u, p, t)
    intidx, nintidx, bemsliptotractotal, els, eta, thetalaw, dc, blockvxglobal, blockvyglobal, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt = p
    @views multmatvec!(vx, vy, els.rotmat[intidx, :, :], u[1:3:end], u[2:3:end])
    @views dtracglobaldt = bemsliptotractotal * [blockvxglobal .- u[1:3:end] blockvyglobal .- u[2:3:end]]'[:]
    @views multmatvec!(dtracxglobaldt, dtracyglobaldt, els.rotmat[intidx, :, :], dtracglobaldt[1:2:end], dtracglobaldt[2:2:end])
    for i in 1:nintidx
        @views dthetadt[i] = thetalaw(abs(vx[i]), abs(u[3:3:end][i]), dc)
        @views dvxdt[i] = 1 / (eta / els.normalstress[intidx[i]] + els.a[intidx[i]] / abs(vx[i])) * (dtracxglobaldt[i] / els.normalstress[intidx[i]] - els.b[intidx[i]] * dthetadt[i] / u[3:3:end][i])
        dvydt[i] = 0
    end
    @views multmatvec!(dudt[1:3:end], dudt[2:3:end], els.rotmatinv[intidx, :, :], dvxdt, dvydt)
    dudt[3:3:end] = dthetadt
    return nothing
end

# function ex_qdstepandsaveflat()
#     # Constants
#     nsteps = 5000
#     nfreesurf = 100
#     nfault = 200
#     printstep = 100
#     amplitude = 1.0
#     outfilename = string(now()) * ".jld2"
#     siay = 365.25 * 24 * 60 * 60
#     tspan = (0, siay * 1000)
#     abstol = 1e-4
#     reltol = 1e-4
#     mu = 3e10
#     nu = 0.25
#     rho = 2700.0
#     eta = mu / (2.0 * sqrt(mu / rho))
#     dc = 0.05
#     blockvelx = 1e-9
#     blockvely = 1e-9

#     @time faultsliptosurfacedisp = inv(partialsconst["trac"]["freesurftopo"]["freesurftopo"]) * partialsconst["trac"]["fault"]["freesurftopo"]
#     faultsliptosurfacedisptofaulttraction = partialsconst["trac"]["freesurftopo"]["fault"] * faultsliptosurfacedisp
#     bemsliptotractotal = partialsconst["trac"]["fault"]["fault"] - faultsliptosurfacedisptofaulttraction
    
#     # CS elements - Euler style stress integration
#     nnodes = 1 * nfault
#     ics = zeros(3 * nnodes)
#     ics[1:3:end] = 1e-3 * blockvelx * ones(nnodes)
#     ics[2:3:end] = 1e-3 * blockvely * ones(nnodes)
#     ics[3:3:end] = 1e8 * ones(nnodes)
#     intidx = idx["fault"] # indices of elements to integrate
#     nintidx = length(intidx)
#     dthetadt = zeros(nintidx)
#     dvxdt = zeros(nintidx)
#     dvydt = zeros(nintidx)
#     vx = zeros(nintidx)
#     vy = zeros(nintidx)
#     dtracxglobaldt = zeros(nintidx)
#     dtracyglobaldt = zeros(nintidx)
#     dtracglobaldt = zeros(2 * nintidx)
#     p = (intidx, nintidx, bemsliptotractotal, els, eta, thetaaginglaw, dc, blockvelx, blockvely, dthetadt, dvxdt, dvydt, vx, vy, dtracxglobaldt, dtracyglobaldt, dtracglobaldt)
#     prob = ODEProblem(derivsconst!, ics, tspan, p)
#     integrator = init(prob, Vern7(), abstol = abstol, reltol = reltol)
#     @time for i in 1:nsteps
#         step!(integrator)
#         if mod(i, printstep) == 0
#             println("step: " * string(i) * " of " * string(nsteps) * ", time: " * string(integrator.sol.t[end] / siay))
#         end    
#     end
#     plotqdtimeseries(integrator.sol, 3, nfault)
    
#     @time @save outfilename integrator.sol els mu nu
#     println("Wrote integration results to:")
#     println(outfilename)
# end
# ex_qdstepandsaveflat()





