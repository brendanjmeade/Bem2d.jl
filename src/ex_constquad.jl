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
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, 4)
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
    # dispconst = zeros(length(xobs), 2)
    # stressconst = zeros(length(xobs), 3)
    # for i in 1:elements.endidx
    #     disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[i],
    #         mu, nu, 0.0, 1.0, elements.xcenter[i], elements.ycenter[i],
    #         elements.rotmat[i, :, :], elements.rotmatinv[i, :, :])
    #     dispconst += disp
    #     stressconst += stress
    # end
    # plotfields(elements, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
    #     dispconst, stressconst, "const slip")

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

    # Calculate all element to element partials
    # TODO: Redo this with arbitrary src and obs indices rather than in order
    # findall(x -> x == "fault01", elements.name)
    srcstart = 1
    srcend = elements.endidx
    obsstart = 1
    obsend = elements.endidx
    partials_disp = zeros(2 * (srcend - srcstart + 1), 2 * (obsend - obsstart + 1))
    partials_stress = zeros(3 * (srcend - srcstart + 1), 2 * (obsend - obsstart + 1))
    partials_trac = zeros(2 * (srcend - srcstart + 1), 2 * (obsend - obsstart + 1))
    for isrc in srcstart:srcend
        for iobs in obsstart:obsend
            pd, ps, pt = partials_constslip(elements, isrc, iobs, mu, nu)
            partials_disp[2 * (isrc - srcstart) + 1:2 * (isrc - srcstart) + 2,
                          2 * (iobs - obsstart) + 1:2 * (iobs - obsstart) + 2] = pd
            partials_stress[3 * (isrc - srcstart) + 1:3 * (isrc - srcstart) + 3,
                            2 * (iobs - obsstart) + 1:2 * (iobs - obsstart) + 2] = ps
            partials_trac[2 * (isrc - srcstart) + 1:2 * (isrc - srcstart) + 2,
                          2 * (iobs - obsstart) + 1:2 * (iobs - obsstart) + 2] = pt
        end
    end

    return nothing
end
ex_constquad()
