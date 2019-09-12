using Revise
using Bem2d

function ex_constquad()
    μ, ν = 3e10, 0.25

    # Create a flat fault
    els = Elements(Int(1e5))
    nels = 20
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, nels)
    els.x1[els.endidx + 1:els.endidx + nels] = x1
    els.y1[els.endidx + 1:els.endidx + nels] = y1
    els.x2[els.endidx + 1:els.endidx + nels] = x2
    els.y2[els.endidx + 1:els.endidx + nels] = y2
    els.name[els.endidx + 1:els.endidx + nels] .= "fault"
    standardize_elements!(els)

    # Set fault slip
    xcenters = els.xcenter[1:els.endidx]
    ycenters = els.ycenter[1:els.endidx]
    xnodes = sort(els.xnodes[1:els.endidx, :][:])
    ynodes = sort(els.ynodes[1:els.endidx, :][:])
    srcidx = findall(x->x == "fault", els.name)

    # # Constant x-slip only
    # constxslip = ones(nels)
    # constyslip = zeros(nels)
    # quadxslip = repeat(constxslip, 1, 3)
    # quadyslip = repeat(constyslip, 1, 3)

    # Linear x-slip only
    slope = 0.001
    constxslip = slope .* xcenters
    constyslip = zeros(nels)
    quadxslip = transpose(reshape(slope .* xnodes, 3, nels))
    quadyslip = zeros(size(quadxslip))

    # Observation coordinates for far-field calculation
    npts = 50; obswidth = 20e3
    xobs, yobs = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Displacements and stresses
    uconst, σconst = constuσ(slip2uσ, xobs, yobs, els, srcidx, constxslip, constyslip, μ, ν)
    uquad, σquad = quaduσ(slip2uσ, xobs, yobs, els, srcidx, quadxslip, quadyslip, μ, ν)

    # Plot
    close("all")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uconst, σconst, "constant slip - constant slip element")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uquad, σquad, "constant slip - quadratic slip element")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uconst - uquad, σconst - σquad, "residuals")
    show()

    return nothing
end
ex_constquad()
