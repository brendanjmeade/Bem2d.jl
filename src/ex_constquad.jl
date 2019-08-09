using Plots
using Revise
using Bem2d
pyplot()

function ex_constquad()
    L = 10e3
    mu = 3e10
    nu = 0.25

    # Create a flat fault
    elements = Elements()
    x1, y1, x2, y2 = discretizedline(-L, -L, L, L, 1)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "fault01"
    end
    standardize_elements!(elements)

    # Observation coordinates for far-field calculation
    npts = 20
    width = 20e3
    xobs = LinRange(-width, width, npts)
    yobs = LinRange(-width, width, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # TODO: Inclined fault
    # TODO: Quadratic element comparision
    # Constant slip element
    dispconst = zeros(length(xobs), 2)
    stressconst = zeros(length(xobs), 3)
    for i in 1:elements.endidx
        disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[i],
            mu, nu, sqrt(2)/2, sqrt(2)/2, elements.xcenter[i], elements.ycenter[i],
            elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
        dispconst += disp
        stressconst += stress
    end
    plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        dispconst, stressconst, "const slip")

    # Constant traction element
    # dispconst = zeros(length(xobs), 2)
    # stressconst = zeros(length(xobs), 3)
    # for i in 1:elements.endidx
    #    disp, stress = dispstress_consttrac(xobs, yobs, elements.halflength[i],
    #                   mu, nu, 1.0, 0.0, elements.xcenter[i], elements.ycenter[i],
    #                   elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
    #    dispconst += disp
    #    stressconst += stress
    # end
    # plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
    #            dispconst, stressconst, "const traction")

    # Calculate constant slip partial derivatives
    # srcidx = findall(x -> x == "fault01", elements.name)
    # obsidx = findall(x -> x == "fault01", elements.name)
    # partials_disp, partials_stress, partials_trac = partials_constslip(elements, srcidx, obsidx, mu, nu)
    return nothing
end
ex_constquad()
