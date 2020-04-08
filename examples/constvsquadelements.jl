using Revise
using Bem2d


"""
    constvsquadelements()

Compare farfield displacements and stresses from constant displacement
and quadratic displacement elements.  The point of this is to highlight 
the fact the quadratic elements can be used to create linear slip elements
that can accurately model a linear displacement gradient with no spurious
singularities.
"""
function constvsquadelements()
    mu = 3e10
    nu = 0.25

    # Create a flat fault
    els = Elements(Int(1e5))
    nels = 20
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, nels)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Convenience dictionary for element names
    idx = getidxdict(els)

    # Set fault slip
    xcenters = els.xcenter[1:els.endidx]
    ycenters = els.ycenter[1:els.endidx]
    xnodes = sort(els.xnodes[1:els.endidx, :][:])
    ynodes = sort(els.ynodes[1:els.endidx, :][:])

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
    npts = 100
    obswidth = 3e3
    xobs, yobs = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Displacements and stresses
    dispconst, stressconst = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], constxslip, constyslip, mu, nu)
    dispquad, stressquad = quaddispstress(slip2dispstress, xobs, yobs, els, idx["fault"], quadxslip, quadyslip, mu, nu)

    # Plot
    close("all")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
               dispconst, stressconst, "constant slip - constant slip element")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
               dispquad, stressquad, "constant slip - quadratic slip element")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
               dispconst - dispquad, stressconst - stressquad, "residuals")
    show()

    return nothing
end
constvsquadelements()
