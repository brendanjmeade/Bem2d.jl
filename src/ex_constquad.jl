using Revise
using Bem2d

function ex_constquad()
    L = 10e3
    μ = 3e10
    ν = 0.25

    # Create a flat fault
    els = Elements()
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nfault)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Observation coordinates for far-field calculation
    npts = 30
    obswidth = 20e3
    xobs, yobs = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Constant slip with constant element
    uconst, σconst = constuσ(slip2uσ, xobs, yobs, els, 1, 1, 0, μ, ν)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uconst, σconst, "constant slip element")

    # Constant slip with quadratic elements
    uquad, σquad = quadslip(xobs, yobs, els, 1, [1 1 1], [0 0 0], μ, ν)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uquad, σquad, "quadratic slip element")

    return nothing
end
ex_constquad()
