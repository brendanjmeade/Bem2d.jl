using Revise
using Statistics
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function discretized_arc(θstart, θend, radius, n_pts)
    # Create geometry of discretized arc
    θrange = collect(LinRange(θstart, θend, n_pts + 1))
    x = @. radius * cos(θrange)
    y = @. radius * sin(θrange)
    x1 = x[1:1:end-1]
    x2 = x[2:1:end]
    y1 = y[1:1:end-1]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end

function calc_brazil(p, r, θ, θ0, R)
    mmax = 1000 # Max number of terms in Hondros series

    #! Analytic stresses in cylindrical coordinates
    σrr = zeros(length(r))
    σθθ = zeros(length(r))
    σrθ = zeros(length(r))
    σrrconstterm = 2.0 * θ0 * -p / deg2rad(180)
    σθθconstterm = 2.0 * θ0 * -p / deg2rad(180)
    leadingterm = 2.0 * -p / deg2rad(180)
    for m in 1:mmax
        σrr += @. (r/R)^(2*m-2) * (1-(1-1/m)*(r/R)^2) * sin(2*m*θ0) * cos(2*m*θ)
        σθθ += @. (r/R)^(2*m-2) * (1-(1+1/m)*(r/R)^2) * sin(2*m*θ0) * cos(2*m*θ)
        σrθ += @. ((r/R)^(2*m) - (r/R)^(2*m-2)) * sin(2*m*θ0) * sin(2*m*θ)
    end
    σrr = @. σrrconstterm .+ leadingterm .* σrr
    σθθ = @. σθθconstterm .- leadingterm .* σθθ
    σrθ = @. leadingterm .* σrθ

    #! Convert analytic cylindrical stresses to Cartesian
    σxx = zeros(length(r))
    σyy = zeros(length(r))
    σxy = zeros(length(r))
    for i in 1:length(r) # Project a single matrix to Cartesian coordinates
        cylindrical_stress_tensor = [σrr[i] σrθ[i] ; σrθ[i] σθθ[i]]
        transformation_matrix = [cos(θ[i]) -sin(θ[i]) ; sin(θ[i]) cos(θ[i])]
        cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transpose(transformation_matrix)
        σxx[i] = cartesian_stress_tensor[1, 1]
        σyy[i] = cartesian_stress_tensor[2, 2]
        σxy[i] = cartesian_stress_tensor[1, 2]
    end
    return σrr, σθθ, σrθ, σxx, σyy, σxy
end

# function calc_brazil(P, x, y, θ0, R)
#     #! Taken from appendix A of:
#     #! https://www.mdpi.com/1996-1073/11/2/304/htm

#     # T is again some weird little scaleing haveing to do with the length of 
#     # the compressed cylinder
#     t = 1
#     r1 = zeros(size(x))
#     r2 = zeros(size(x))
#     theta1 = zeros(size(x))
#     theta2 = zeros(size(x))

#     #! Analytic stresses in Cartesian coordinates
#     Sxx = zeros(length(x))
#     Syy = zeros(length(x))
#     Sxy = zeros(length(x))
#     for i in 1:length(x) # Project a single matrix to Cartesian coordinates
#         Sxx[i] = cartesian_stress_tensor[1, 1]
#         Syy[i] = cartesian_stress_tensor[2, 2]
#         Sxy[i] = cartesian_stress_tensor[1, 2]
#     end
#     return Sxx, Syy, Sxy
# end


function ex_cylinder_profile()
    PyPlot.close("all")
    mu = 3e10
    nu = 0.25

    #! Try Crouch and Starfield line style plot
    # Figure 4.15 in CS1983 doesn't make a lot of sense to me because
    # it shows a sign change in sigmaxxdivp which is not seen the data given
    # in table 4.1 
    xdivr = [0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
    Sxxdivp = [0.0398 0.0382 0.0339 0.0278 0.0209 0.0144 0.0089 0.0047 0.0019 0.0004 0.000]
    Syydivp = -[0.1198 0.1167 0.1078 0.0946 0.0789 0.0624 0.0465 0.0321 0.0195 0.0089 0.000]

    # Solution from Hondros (1959) as summarized by Wei and Chau 2013
    p = -1.0e5 # Applied radial pressure over arc
    arclen = deg2rad(45.0) # Arc length over which pressure is applied
    R = 57.296 # Radius of disc
    y = @. R * collect(0.0:0.1:1.0)
    x = zeros(size(y))
    r = @. sqrt(x^2 + y^2)
    theta = @. atan(y, x)

    #! Analytic Brazil test
    _, _, _, Sxxhondros, Syyhondros, Sxyhondros = calc_brazil(p, r, theta, arclen, R)

    #! BEM solution
    els = Bem2d.Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, 360)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "circle"
    end
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)
    partialsconst = Bem2d.initpartials(els)

    #! Apply normal tractions everywhere and convert from radial to Cartesian
    xtrac = zeros(els.endidx)
    ytrac = zeros(els.endidx)
    thetaels = @. atan(els.ycenter[1:1:els.endidx], els.xcenter[1:1:els.endidx])
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtrac[i], ytrac[i] = els.rotmat[i, :, :] * normalTractions
    end

    #! Zero out the tractions on the area without contact
    deleteidx = findall(x -> (x>arclen && x<deg2rad(180)-arclen), thetaels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0
    deleteidx = findall(x -> (x<-arclen && x>-deg2rad(180)+arclen), thetaels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0

    #! Scale tractions by element lengths
    xtracscaled = xtrac ./ els.length[1:1:els.endidx]
    ytracscaled = ytrac ./ els.length[1:1:els.endidx]

    # Given the other tractions calcuate the induced displacements on the boudaries
    T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["circle"], idx["circle"], mu, nu)
    U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["circle"], idx["circle"], mu, nu)
    # dispall = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * Bem2d.interleave(xtrac, ytrac)
    dispall = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * Bem2d.interleave(xtracscaled, ytracscaled)

    #! Streses from tractions and displacements
    _, Strac = Bem2d.constdispstress(trac2dispstress, x, y, els, idx["circle"], xtracscaled, ytracscaled, mu, nu)
    _, Sdisp = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["circle"], dispall[1:2:end], dispall[2:2:end], mu, nu)
    Sbem = @. Strac + Sdisp

    # Plot profiles
    fontsize = 24
    markersize = 20
    figure(figsize=(15,15))

    # xx component of stresses
    ax = subplot(2, 1, 1)
    ax.tick_params("both", labelsize=fontsize)
    # for i in 1:length(xdivr)
    #     if i == 1
    #         plot(xdivr[i], Sxxdivp[i], "rs", markersize=markersize, label="CS table 4.1")
    #     else
    #         plot(xdivr[i], Sxxdivp[i], "rs", markersize=markersize)
    #     end
    # end
    plot(y[1:end-1]./R, -Sxxhondros[1:end-1] ./ p, "g^", markersize=markersize, label="Hondros")
    plot(y[1:end-1]./R, -Sbem[1:end-1, 1] ./ p, "bs", markersize=markersize, label="BEM")
    legend(fontsize=fontsize)
    xlabel(L"x/R", fontsize=fontsize)
    ylabel(L"\sigma_{xx}/p", fontsize=fontsize)

    #! yy component of stresses
    ax = subplot(2, 1, 2)
    ax.tick_params("both", labelsize=fontsize)
    # for i in 1:length(xdivr)
    #     if i == 1
    #         plot(xdivr[i], Syydivp[i], "rs", markersize=markersize, label="CS table 4.1")
    #     else
    #         plot(xdivr[i], Syydivp[i], "rs", markersize=markersize)
    #     end
    # end
    plot(y[1:end-1]./R, -Syyhondros[1:end-1] ./ p, "g^", markersize=markersize, label="Hondros")
    plot(y[1:end-1]./R, -Sbem[1:end-1, 2] ./ p, "bs", markersize=markersize, label="BEM")
    legend(fontsize=fontsize)
    xlabel(L"x/R", fontsize=fontsize)
    ylabel(L"\sigma_{yy}/p", fontsize=fontsize)
    show()
end
ex_cylinder_profile()
