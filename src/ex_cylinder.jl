import Bem2d
import Infiltrator

function ex_cylinder()

    # Observation coordinates for far-field calculation
    npts = 20
    obswidth = 2e3
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Convert Cartesian to cylindrical coordinates
    r = @. sqrt(x^2 + y^2)
    θ = @. atan(y, x)
            
    # Solution from Hondros (1959) as summarized by Wei and Chau 2013
    σrr = zeros(length(x))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))
    p = 1.0e5 # Applied radial pressure over arc
    θ0 = 90.0 # Arc length over which pressure is applied
    R = 1.0e3 # Radius of disc
    mmax = 100 # Max number of terms in Hondros series

    for m in 1:mmax

    end
    
    # Convert cylindrical stresses to Cartesian
    # http://solidmechanics.org/text/AppendixD/AppendixD.htm

end
ex_cylinder()
