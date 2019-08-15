using Revise
using Bem2d

function ex_thrusttopo()
    mu = 30e9
    nu = 0.25
    elements = Elements()

    # Observation points for internal evaluation and visualization
    npts= 100
    xobs = LinRange(-10e3, 10e3, npts)
    yobs = LinRange(-5e3, 5e3, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # Topographic free surface
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, 20)
    y1 = -1e3 * @.atan(x1 / 1e3)
    y2 = -1e3 * @.atan(x2 / 1e3)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "freesurface"
    end
    standardize_elements!(elements)

    # Curved fault
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, 10)
    y1 = 3e3 * @.atan(x1 / 1e3)
    y2 = 3e3 * @.atan(x2 / 1e3)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "fault"
    end
    standardize_elements!(elements)

    # Partial derivatves
    srcidx = findall(x -> x == "fault", elements.name)
    obsidx = findall(x -> x == "freesurface", elements.name)
    d1, s1, t1 = partials_constslip(elements, srcidx, obsidx, mu, nu)
    srcidx = findall(x -> x == "freesurface", elements.name)
    obsidx = findall(x -> x == "freesurface", elements.name)
    d2, s2, t2 = partials_constslip(elements, srcidx, obsidx, mu, nu)

    # Remove and separate BCs, local to global transform here?
    nfaultelements = length(findall(x -> x == "fault", elements.name))
    faultslip = zeros(2 * nfaultelements)
    faultslip[1:2:end] .= 1.0 # Global coordinate system

    # Solve the BEM problem
    disp_freesurface = inv(t2) * t1 * faultslip

    # Fault in full space
    dispfault = zeros(length(xobs), 2)
    stressfault = zeros(length(xobs), 3)
    faultidx = findall(x -> x == "fault", elements.name)
    for i in 1:length(faultidx)
        disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[faultidx[i]],
            mu, nu, 1, 0, elements.xcenter[faultidx[i]], elements.ycenter[faultidx[i]],
            elements.rotmat[faultidx[i], :, :], elements.rotmatinv[faultidx[i], :, :])
        dispfault += disp
        stressfault += stress
    end
    # plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
    #     dispfault, stressfault, "fault")

    # Free surface in full space
    dispfreesurface::Array{Float64} = zeros(length(xobs), 2)
    stressfreesurface::Array{Float64} = zeros(length(xobs), 3)
    freesurfaceidx::Array{Int64} = findall(x -> x == "freesurface", elements.name)
    for i in 1:length(freesurfaceidx)
        disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[freesurfaceidx[i]],
            mu, nu, disp_freesurface[1:2:end][i], disp_freesurface[2:2:end][i],
            elements.xcenter[freesurfaceidx[i]], elements.ycenter[freesurfaceidx[i]],
            elements.rotmat[freesurfaceidx[i], :, :], elements.rotmatinv[freesurfaceidx[i], :, :])
        dispfreesurface += disp
        stressfreesurface += stress
    end
    # plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
    #     dispfreesurface, stressfreesurface, "free surface")

    # Plot fault + free surface
    plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        dispfault + dispfreesurface, stressfault + stressfreesurface, "total")
end
ex_thrusttopo()
