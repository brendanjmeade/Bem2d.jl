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
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, 3)
    for i in 1:length(x1)
        elements.x1[i + elements.lastidx] = x1[i]
        elements.y1[i + elements.lastidx] = y1[i]
        elements.x2[i + elements.lastidx] = x2[i]
        elements.y2[i + elements.lastidx] = y2[i]
        elements.uxconst[i + elements.lastidx] = 1.0
        elements.uyconst[i + elements.lastidx] = 1.0
        elements.uxquad[i + elements.lastidx, :] = [1.0 1.0 1.0]
        elements.uyquad[i + elements.lastidx, :] = [1.0 1.0 1.0]
        elements.name[i + elements.lastidx] = "fault01"
    end
    standardize_elements!(elements)

    # Observation coordinates for far-field calculation
    npts = 100
    width = 20e3
    xobs = range(-width, stop=width, length=npts)
    yobs = range(-width, stop=width, length=npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # A simple forward model for the volume
    dispconst = zeros(length(xobs), 2)
    stressconst = zeros(length(xobs), 3)
    for i in 1:elements.lastidx
        disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[i],
                       mu, nu, 1.0, 0.0, elements.xcenter[i], elements.ycenter[i],
                       elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
        dispconst += disp
        stressconst += stress
    end

    plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
               dispconst, stressconst, "constant elements (slip)")

    return nothing
end

ex_constquad()
