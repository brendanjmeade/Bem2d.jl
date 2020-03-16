using Revise
using Infiltrator
using PyPlot

function stressfunction(slip, mu, x, y, d1, d2)
    stressd1 = @. slip*mu/(2*pi) * ((y+d1)/(x^2+(y+d1)^2) - (y+d1)/(x^2+(y-d1)^2))
    stressd2 = @. slip*mu/(2*pi) * ((y+d2)/(x^2+(y+d2)^2) - (y+d2)/(x^2+(y-d2)^2))
    stress = stressd1 .- stressd2
    return stress
end

function stressforplotting(stress)
    stress = @. sign(stress) * abs(stress)^(1.0/3.0)
    return stress
end

function interactionmatrix(nels, dtop, dbot, mu)
    mat = zeros(nels, nels) # to store results
    dvec = collect(LinRange(dtop, dbot, nels+1))
    dtopvec = dvec[1:1:end-1]
    dbotvec = dvec[2:1:end]
    ymid = @. -(dtopvec + dbotvec) / 2 
    
    # Loop over all possible combinations
    for i = 1:nels
        for j = 1:nels
            mat[i, j] = stressfunction(1, mu, 0, ymid[i], dtopvec[j], dbotvec[j])
        end
    end
    return mat
end


function strikeslipstress()
    PyPlot.close("all")
    mu = 3e10
    d1 = 20e3 # meters
    d2 = 0e3 # meters
    maxdepth = -50e3 # meters
    npts = 1000
    x = @. 0 * ones(npts)
    y = collect(LinRange(0, maxdepth, npts))
    stressco = stressfunction(1, mu, x, y, d1, d2)
    stressinter = stressfunction(1, mu, x, y, 1000e3, 25e3)

    # Calculate full element to element interaction interaction
    mat = interactionmatrix(10, 0e3, 20e3, mu)

    PyPlot.matshow(mat) # Plot interaction matrix
    yplot = y / 1e3 # convert y to km for plotting convenience
    fontsize = 24
    PyPlot.figure(figsize=(10,15))
    PyPlot.plot(stressforplotting(stressco), yplot, "-r")
    PyPlot.plot(stressforplotting(stressinter), yplot, "-b")
    PyPlot.ylabel("y (m)", fontsize=fontsize)
    PyPlot.xlabel("shear (Pa)", fontsize=fontsize)
    gca().set_ylim([minimum(yplot), 0])
    gca().fontsize= 30

    PyPlot.show()

end
strikeslipstress()




