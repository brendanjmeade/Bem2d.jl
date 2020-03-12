using Revise
using Test
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function pol2cart(theta, radius)
    x = @. radius * cosd(theta)
    y = @. radius * sind(theta)
    return x, y
end

function ex_single()
    mu = 3e10
    nu = 0.25
    atol = 1e-20 # absolute tolerance for element wise testing

    # Forcing (traction and displacement)
    appliedx = [1.0 ; 0.0]
    appliedy = [0.0 ; 1.0]

    # Angle to x, y
    thetavec = collect(0:45:90)
    theta = 90 * rand(1)[1]
    theta = 45
    x, y = pol2cart(theta, 0.5)

    # Elements
    els = Bem2d.Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-x, -y, x, y, 1)
    els.x1[1] = x1[1]
    els.y1[1] = y1[1]
    els.x2[1] = x2[1]
    els.y2[1] = y2[1]
    els.name[1] = "line"
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)
    
    # All the matrices for the single element
    T, S, H = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["line"], idx["line"], mu, nu)
    U, D, A = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["line"], idx["line"], mu, nu)
    
    # Compare forward model at element centroid with partials
    # These should agree...and if they do that's a problem!!!
    disptrac, stresstrac = Bem2d.constdispstress(trac2dispstress, 0, 0, els, idx["line"], appliedx[1], appliedx[2], mu, nu)
    @show disptrac
    
    println("\n *** theta = " * string(theta) * " (degrees) ***")
    println("\n *** T ***")
    display(T)
    println("\n *** U ***")
    display(U)
    @show(U)
    # Tractions -> displacements -> tractions (x-forcing)
    dispinducedx = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * appliedx
    tracrecoveredx = inv(U) * (T + 0.5 * LinearAlgebra.I(size(T)[1])) * dispinducedx
    println("\n *** Tractions -> displacements -> tractions (x-forcing) ***")
    display(@test all(isapprox.(appliedx, tracrecoveredx, atol=atol)))

    # Tractions -> displacements -> tractions (y-forcing)
    dispinducedy = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * appliedy
    tracrecoveredy = inv(U) * (T + 0.5 * LinearAlgebra.I(size(T)[1])) * dispinducedy
    println("\n *** Tractions -> displacements -> tractions (y-forcing) ***")
    display(@test all(isapprox.(appliedy, tracrecoveredy, atol=atol)))

    # Displacements -> tractions -> displacements (x-forcing)
    tracinducedx = inv(U) * (T + 0.5 * LinearAlgebra.I(size(T)[1]))* appliedx 
    disprecoveredx =  (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * tracinducedx
    println("\n *** Displacements -> tractions -> displacements (x-forcing) ***")
    display(@test all(isapprox.(appliedx, disprecoveredx, atol=atol)))

    # Displacements -> tractions -> displacements (x-forcing)
    tracinducedy = inv(U) * (T + 0.5 * LinearAlgebra.I(size(T)[1]))* appliedy 
    disprecoveredy =  (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * tracinducedy
    println("\n *** Displacements -> tractions -> displacements (y-forcing) ***")
    # println("??? " * string(:appliedy) * " == " * string(:disprecoveredy) * " ???")
    display(@test all(isapprox.(appliedy, disprecoveredy, atol=atol)))
end
ex_single()
