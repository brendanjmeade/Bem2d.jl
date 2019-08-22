using Revise
using Bem2d
# using Bem2dQuadKernels

function ex_constquad()
    L = 10e3
    μ = 3e10
    ν = 0.25

    # Create a flat fault
    els = Elements()
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-L, -L, L, L, nfault)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Observation coordinates for far-field calculation
    npts = 20
    width = 20e3
    xobs = LinRange(-width, width, npts)
    yobs = LinRange(-width, width, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # Constant slip element
    uconst = zeros(length(xobs), 2)
    σconst = zeros(length(xobs), 3)
    for i in 1:els.endidx
        u, σ = constslip(xobs, yobs, els.halflength[i],
            μ, ν, sqrt(2) / 2, sqrt(2) / 2, els.xcenter[i], els.ycenter[i],
            els.rotmat[i, :, :], els.rotmatinv[i, :, :])
        uconst += u
        σconst += σ
    end
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uconst, σconst, "constant slip element")

    # TODO: Should I push the loop over elements into quadslip/constslip?
    #   It would behave the same for a single element
    # Quadratic elements
    uquad = zeros(length(xobs), 2)
    σquad = zeros(length(xobs), 3)
    for i in 1:els.endidx
        u, σ = quadslip(xobs, yobs, els.halflength[i], μ, ν,
            [sqrt(2) / 2 sqrt(2) / 2 sqrt(2) / 2], [sqrt(2) / 2 sqrt(2) / 2 sqrt(2) / 2],
            els.xcenter[i], els.ycenter[i],
            els.rotmat[i, :, :], els.rotmatinv[i, :, :])

        # u, σ = quadslip(xobs, yobs, els, idx, μ, ν, [1 1 1], [0 0 0])
        uquad += u
        σquad += σ
    end
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uquad, σquad, "quadratic slip element")

    return nothing
end
ex_constquad()
