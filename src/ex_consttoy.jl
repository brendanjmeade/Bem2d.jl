using Revise
using LinearAlgebra
using PyPlot
using Bem2d


function ex_consttoy()
    PyPlot.close("all")
    mu = 3e10
    nu = 0.25

    # Create a flat fault
    els = Elements(Int(1e5))
    nels = 1
    x1, y1, x2, y2 = discretizedline(-5e3, 0, 5e3, 0, nels)
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
    npts = 200
    obswidth = 20e3
    xobs, yobs = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)


    # Volume evalattion of Crouch and Starfield (1983) kernels
    dispconstslip, stressconstslip = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    dispconsttrac, stressconsttrac = constdispstress(trac2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts), dispconstslip, stressconstslip, "T, S (CS displacement kernels)")
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts), dispconsttrac, stressconsttrac, "U, D (CS stress kernels)")
    show()

    # Constant traction element
    # Generate partials
    T, _, _ = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)
    U, _, _ = partialsconstdispstress(trac2dispstress, els, idx["fault"], idx["fault"], mu, nu)

    # Solve BEM problem for slip on fault resulting from unit traction
    u = (inv(T + 0.5 * I(size(T)[1]))) * U * [xdrive; ydrive]

    # Forward model for volume
    uconstslip, sconstslip = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], u[1], u[2], mu, nu)
    uconsttrac, sconsttrac = constdispstress(trac2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
               uconsttrac + uconstslip, sconsttrac + sconstslip,
               "constant traction element")

    # Constant displacement element
    # Generate partials - same as above
    # Solve BEM problem for slip on fault resulting from unit traction
    t = inv(U) * (T + 0.5 * I(size(T)[1])) * [xdrive; ydrive]

    # Forward model for volume
    # uconstslip, sconstslip = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], t[1], t[2], mu, nu)
    # uconsttrac, sconsttrac = constdispstress(trac2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    uconstslip, sconstslip = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    uconsttrac, sconsttrac = constdispstress(trac2dispstress, xobs, yobs, els, idx["fault"], t[1], t[2], mu, nu)

    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts),
               uconsttrac + uconstslip, sconsttrac + sconstslip,
               "constant displacement element")


    ### Displacement discontinuity element
    dispconstslip, stressconstslip = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], xdrive, ydrive, mu, nu)
    plotfields(els, reshape(xobs, npts, npts), reshape(yobs, npts, npts), dispconstslip, stressconstslip, "Displacement discontinuity element")
    PyPlot.show()

    return nothing
end
ex_consttoy()
