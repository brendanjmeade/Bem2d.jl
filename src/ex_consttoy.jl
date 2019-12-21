using Revise
using LinearAlgebra
using Bem2d

function ex_consttoy()
    mu = 3e10
    nu = 0.25

    # Create a flat fault
    els = Elements(Int(1e5))
    nels = 1
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
    par = initpartials(els)
    
    # Set driving traction / slip
    xdrive = ones(nels)
    ydrive = zeros(nels)

    # Observation coordinates for far-field calculation
    npts = 100
    obswidth = 20e3
    xobs, yobs = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    ### Show Crouch and Starfield (1983) kernels
    dispconstslip, stressconstslip = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    dispconsttrac, stressconsttrac = constdispstress(trac2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts), dispconstslip, stressconstslip, "CS slip kernel only")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts), dispconsttrac, stressconsttrac, "CS stress/traction kernel only")
    show()

    ### Constant traction element
    # Generate partials
    T, _, _ = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    U, _, _ = partialsconstdispstress(trac2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    
    # Solve BEM problem for slip on fault resulting from unit traction
    u = (I(size(T)[1]) - T * inv(T - I(size(T)[1]))) * U * [xdrive; ydrive]
    @show u
    
    # TODO: Forward model for volume

    ### Constant displacement element
    # TODO: Generate partials
    # TODO: Solve BEM problem
    # TODO: Forward model for volume

    ### Displacement discontinuity element
    # TODO: Generate partials
    # TODO: Solve BEM problem
    # TODO: Forward model for volume
    
    return nothing
end
ex_consttoy()
