using Revise
using Bem2d

function ex_constquad()
    μ, ν = 3e10, 0.25

    # Create a flat fault
    els = Elements(Int(1e5))
    nels = 1
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, nels)
    els.x1[els.endidx + 1:els.endidx + nels] = x1
    els.y1[els.endidx + 1:els.endidx + nels] = y1
    els.x2[els.endidx + 1:els.endidx + nels] = x2
    els.y2[els.endidx + 1:els.endidx + nels] = y2
    els.name[els.endidx + 1:els.endidx + nels] .= "fault"
    standardize_elements!(els)

    # Observation coordinates for far-field calculation
    npts = 50; obswidth = 20e3
    xobs, yobs = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Constant slip with constant element
    uconst, σconst = constuσ(slip2uσ, xobs, yobs, els, 1, 1, 0, μ, ν)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uconst, σconst, "constant slip - constant slip element")

    # Constant slip with quadratic elements
    uquad, σquad = quaduσ(slip2uσ, xobs, yobs, els, 1, [1 1 1], [0 0 0], μ, ν)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        uquad, σquad, "constant slip - quadratic slip element")

    return nothing
end
ex_constquad()
