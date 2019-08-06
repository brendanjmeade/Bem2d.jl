using Plots
using Revise
using Bem2d
pyplot()

function ex_constquad()
    L = 5e3
    mu = 3e10
    nu = 0.25

    # Create a flat fault
    elements = Elements()
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, 1)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.uxconst[i + elements.endidx] = 1.0
        elements.uyconst[i + elements.endidx] = 1.0
        elements.uxquad[i + elements.endidx, :] = [1.0 1.0 1.0]
        elements.uyquad[i + elements.endidx, :] = [1.0 1.0 1.0]
        elements.name[i + elements.endidx] = "fault01"
    end
    standardize_elements!(elements)
    println(elements.endidx)
    println(fieldnames(Elements))

    # Observation coordinates for far-field calculation
    npts = 50
    width = 20e3
    xobs = range(-width, stop=width, length=npts)
    yobs = range(-width, stop=width, length=npts)
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
                       mu, nu, 1.0, 0.0, elements.xcenter[i], elements.ycenter[i],
                       elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
        dispconst += disp
        stressconst += stress
    end
    # plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
    #            dispconst, stressconst, "const slip")

    # Constant traction element
    dispconst = zeros(length(xobs), 2)
    stressconst = zeros(length(xobs), 3)
    for i in 1:elements.endidx
       disp, stress = dispstress_consttrac(xobs, yobs, elements.halflength[i],
                      mu, nu, 1.0, 0.0, elements.xcenter[i], elements.ycenter[i],
                      elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
       dispconst += disp
       stressconst += stress
    end
    # plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
    #            dispconst, stressconst, "const traction")

    xobs = elements.xcenter
    yobs = elements.ycenter
    normalvector = [elements.xnormal[1] ; elements.ynormal[1]]
    pd, ps, pt = partials_constslip(elements, 1, xobs, yobs, normalvector, mu, nu)
    println(pd)
    println(ps)
    println(pt)

    return nothing
end
ex_constquad()
