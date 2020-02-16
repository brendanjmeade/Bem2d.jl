using Revise
using LaTeXStrings
using PyPlot
using LinearAlgebra
using Infiltrator
using Bem2d

function discretized_arc(θstart, θend, radius, n_pts)
    # Create geometry of discretized arc
    θrange = collect(LinRange(θstart, θend, n_pts + 1))
    x = @. radius * cosd(θrange)
    y = @. radius * sind(θrange)
    x1 = x[1:1:end-1]
    x2 = x[2:1:end]
    y1 = y[1:1:end-1]
    y2 = y[2:1:end]
    return x1, y1, x2, y2
end

function circle_subplot(x, y, mat, npts, R, θ0, title_string)
    fontsize = 20
    contour_levels = 20
    contour_color = "white"
    contour_linewidth = 0.5
    color_scale = 1
    # mat = @. mat / 1e6 # Convert from Pascals to mega Pascals

    PyPlot.contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = PyPlot.colorbar(fraction=0.020, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.set_label(label=title_string * " (Pa)", fontsize=fontsize)
    PyPlot.contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)

    # Draw enditre circle
    x1, y1, x2, y2 = discretized_arc(-180, 180, R, 360)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-k", linewidth=10)
    end

    # Draw right-hand side of applied compression arc
    x1, y1, x2, y2 = discretized_arc(-θ0, θ0, R, 50)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-r", linewidth=2)
    end

    # Draw left-hand side of applied compression arc
    x1, y1, x2, y2 = discretized_arc(-θ0+180, θ0+180, R, 50)
    for i in 1:length(x1) # Draw arc where compression is being applied
        PyPlot.plot([x1[i], x2[i]], [y1[i], y2[i]], "-r", linewidth=2)
    end

    PyPlot.xlabel(L"x \; (m)", fontsize=fontsize)
    PyPlot.ylabel(L"y \; (m)", fontsize=fontsize)
    # PyPlot.xlim([-1100, 1100])
    # PyPlot.ylim([-1100, 1100])
    # PyPlot.xticks([-1000, 0, 1000])
    # PyPlot.yticks([-1000, 0, 1000])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)
end

function normalizenan(vec)
    vec = vec .- minimum(filter(!isnan, vec))
    vec = vec ./ maximum(filter(!isnan, vec))
    return vec
end

function ex_cylinder()
    mu = 3e10
    nu = 0.25

    # Solution from Hondros (1959) as summarized by Wei and Chau 2013
    p = -1.0e5 # Applied radial pressure over arc
    θ0 = 10.0 # Arc length over which pressure is applied
    R = 1.0e3 # Radius of disc
    R = 57.296 # Radius of disc
    mmax = 1000 # Max number of terms in Hondros series
    npts = 50
    x, y = Bem2d.obsgrid(-R, -R, R, R, npts)
    r = @. sqrt(x^2 + y^2)
    θ = @. rad2deg(atan(y, x))

    σrr = zeros(length(x))
    σθθ = zeros(length(x))
    σrθ = zeros(length(x))
    σrrconstterm = 2.0 * θ0 * -p / 180
    σθθconstterm = 2.0 * θ0 * -p / 180
    leadingterm = 2.0 * -p / 180
    for m in 1:mmax
        σrr += @. (r/R)^(2*m-2) * (1-(1-1/m)*(r/R)^2) * sind(2*m*θ0) * cosd(2*m*θ)
        σθθ += @. (r/R)^(2*m-2) * (1-(1+1/m)*(r/R)^2) * sind(2*m*θ0) * cosd(2*m*θ)
        σrθ += @. ((r/R)^(2*m) - (r/R)^(2*m-2)) * sind(2*m*θ0) * sind(2*m*θ)
    end
    σrr = @. σrrconstterm + leadingterm * σrr
    σθθ = @. σθθconstterm - leadingterm * σθθ
    σrθ = @. leadingterm * σrθ

    # Convert cylindrical stresses to Cartesian
    σxx = zeros(length(x))
    σyy = zeros(length(x))
    σxy = zeros(length(x))
    for i in 1:length(x) # Project a single matrix to Cartesian coordinates
        cylindrical_stress_tensor = [σrr[i] σrθ[i] ; σrθ[i] σθθ[i]]
        transformation_matrix = [cosd(θ[i]) -sind(θ[i]) ; sind(θ[i]) cosd(θ[i])]
        cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transpose(transformation_matrix)
        σxx[i] = cartesian_stress_tensor[1, 1]
        σyy[i] = cartesian_stress_tensor[2, 2]
        σxy[i] = cartesian_stress_tensor[1, 2]
    end

    # The solution due to line sources at each end of cylinder
    # From: Estimation of the tensile elastic modulus using Brazilian disc by applying
    # diametrically opposed concentrated loads, Ye Jianhong a,, F.Q. Wu , J.Z. Sun
    # Equation 6
    # Sort skeptical about this one because I can't get their (1959) equations to work (Equations 1)
    l = 1.0 # This is the "length" of the cylinder...I have no idea what this should be
    # σxxline = @. 2*p / (pi*l) * (((R-y)*x^2)/((R-y)^2 + x^2)^2 + ((R+y)*x^2)/((R+y)^2 + x^2)^2 - (1/(2*R)))
    # σyyline = @. 2*p / (pi*l) * (((R-y)^3)/((R-y)^2 + x^2)^2 + ((R+y)^3)/((R+y)^2 + x^2)^2 - (1/(2*R)))
    # σxyline = @. 2*p / (pi*l) * (((R-y)^2*x)/((R-y)^2 + x^2)^2 + ((R+y)^2*x)/((R+y)^2 + x^2)^2 - (1/(2*R)))

    # Just swapped (x,y) to (y,x) to make this with line sources at y = 0 rather than x = 0.
    σxxline = @. 2*p / (pi*l) * (((R-x)*y^2)/((R-x)^2 + y^2)^2 + ((R+x)*y^2)/((R+x)^2 + y^2)^2 - (1/(2*R)))
    σyyline = @. 2*p / (pi*l) * (((R-x)^3)/((R-x)^2 + y^2)^2 + ((R+x)^3)/((R+x)^2 + y^2)^2 - (1/(2*R)))
    σxyline = @. 2*p / (pi*l) * (((R-x)^2*y)/((R-x)^2 + y^2)^2 + ((R+x)^2*y)/((R+x)^2 + y^2)^2 - (1/(2*R)))

    # Start of BEM solution
    els = Bem2d.Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(-180, 180, R, 360)
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

    # Apply normal tractions everywhere
    # Convert from radial to Cartesian
    xtrac = zeros(els.endidx)
    ytrac = zeros(els.endidx)
    θels = @. rad2deg(atan(els.ycenter[1:1:els.endidx], els.xcenter[1:1:els.endidx]))
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        normalTractions = [0; p] # Pressure in fault normal component only.
        xtrac[i], ytrac[i] = els.rotmat[i, :, :] * normalTractions
    end

    # Zero out the tractions on the area without contact
    # What are the boundary conditions on the regions without contact?
    deleteidx = findall(x -> (x>θ0 && x<180-θ0), θels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0
    deleteidx = findall(x -> (x<-θ0 && x>-180+θ0), θels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0

    # Given the traction boundary conditions calcuate the induced displacements on each element
    xdisp = zeros(els.endidx)
    ydisp = zeros(els.endidx)
    for i in 1:els.endidx # Calcuate the x and y components of the tractions
        T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["circle"][i], idx["circle"][i], mu, nu)
        U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["circle"][i], idx["circle"][i], mu, nu)
        # xdisp[i], ydisp[i] = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * [xtrac[i]; ytrac[i]]
    end
    dispall_isolated = Bem2d.interleave(xdisp, ydisp)

    # Given the other tractions calcuate the induced displacements on the boudaries
    interleavedtracs = Bem2d.interleave(xtrac, ytrac)
    T, _, _ = Bem2d.partialsconstdispstress(slip2dispstress, els, idx["circle"], idx["circle"], mu, nu)
    U, _, _ = Bem2d.partialsconstdispstress(trac2dispstress, els, idx["circle"], idx["circle"], mu, nu)
    dispall = (inv(T + 0.5 * LinearAlgebra.I(size(T)[1]))) * U * Bem2d.interleave(xtrac, ytrac)
    dispall_Uonly = U * Bem2d.interleave(xtrac, ytrac)

    # Streses from tractions
    _, stresstrac = Bem2d.constdispstress(trac2dispstress, x, y, els, idx["circle"], xtrac, ytrac, mu, nu)

    # Stresses from traction induced displacements
    _, stressdisp = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["circle"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    # Plot tractions
    # Plot the tractions and induced displacements at each element centroid
    # Should the induced displacements be those from *all* tractions? or just a single element
    fontsize = 20
    PyPlot.close("all")
    PyPlot.figure(figsize=(20, 10))
    PyPlot.subplot(1, 2, 1)
    for i in 1:els.endidx
        PyPlot.plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k")
    end
    PyPlot.quiver(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], xtrac, ytrac)
    PyPlot.xlabel(L"x \; (m)", fontsize=fontsize)
    PyPlot.ylabel(L"y \; (m)", fontsize=fontsize)
    PyPlot.title("applied tractions", fontsize=fontsize)
    # PyPlot.xlim([-1100, 1100])
    # PyPlot.ylim([-1100, 1100])
    # PyPlot.xticks([-1000, 0, 1000])
    # PyPlot.yticks([-1000, 0, 1000])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)

    PyPlot.subplot(1, 2, 2)
    for i in 1:els.endidx
        PyPlot.plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k")
    end
    PyPlot.quiver(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], xdisp, ydisp, color="blue")
    PyPlot.quiver(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], dispall[1:2:end], dispall[2:2:end], color="green")
    PyPlot.xlabel(L"x \; (m)", fontsize=fontsize)
    PyPlot.ylabel(L"y \; (m)", fontsize=fontsize)
    PyPlot.title("induced displacements", fontsize=fontsize)
    # PyPlot.xlim([-1100, 1100])
    # PyPlot.ylim([-1100, 1100])
    # PyPlot.xticks([-1000, 0, 1000])
    # PyPlot.yticks([-1000, 0, 1000])
    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().tick_params(labelsize=fontsize)

    # Try setting a few values to NaN and see if we can isolate the circle
    to_nan_idx = findall(x -> x > 1.0 * R, r)
    σrr[to_nan_idx] .= NaN
    σθθ[to_nan_idx] .= NaN
    σrθ[to_nan_idx] .= NaN
    σxx[to_nan_idx] .= NaN
    σyy[to_nan_idx] .= NaN
    σxy[to_nan_idx] .= NaN
    σxxline[to_nan_idx] .= NaN
    σyyline[to_nan_idx] .= NaN
    σxyline[to_nan_idx] .= NaN
    stresstrac[to_nan_idx, 1] .= NaN
    stresstrac[to_nan_idx, 2] .= NaN
    stresstrac[to_nan_idx, 3] .= NaN
    stressdisp[to_nan_idx, 1] .= NaN
    stressdisp[to_nan_idx, 2] .= NaN
    stressdisp[to_nan_idx, 3] .= NaN

    # Try normalizing and sign flips
    σxx = -σxx
    σyy = -σyy
    σxy = -σxy
    σxx = normalizenan(σxx)
    σyy = normalizenan(σyy)
    σxy = normalizenan(σxy)
    stresstrac[:, 1] = normalizenan(stresstrac[:, 1])
    stresstrac[:, 2] = normalizenan(stresstrac[:, 2])
    stresstrac[:, 3] = normalizenan(stresstrac[:, 3])

    # Plot contours with PyPlot
    fontsize = 24
    contour_levels = 10
    contour_color = "black"
    contour_linewidth = 0.5

    # Analytic solutions
    PyPlot.figure(figsize=(20,20))
    PyPlot.subplot(3, 3, 1)
    circle_subplot(x, y, σrr, npts, R, θ0, L"\sigma_{rr}")
    PyPlot.subplot(3, 3, 2)
    circle_subplot(x, y, σθθ, npts, R, θ0, L"\sigma_{\theta\theta}")
    PyPlot.subplot(3, 3, 3)
    circle_subplot(x, y, σrθ, npts, R, θ0, L"\sigma_{r\theta}")
    PyPlot.subplot(3, 3, 4)
    circle_subplot(x, y, σxx, npts, R, θ0, L"\sigma_{xx}")
    PyPlot.subplot(3, 3, 5)
    circle_subplot(x, y, σyy, npts, R, θ0, L"\sigma_{yy}")
    PyPlot.subplot(3, 3, 6)
    circle_subplot(x, y, σxy, npts, R, θ0, L"\sigma_{xy}")
    PyPlot.subplot(3, 3, 7)
    circle_subplot(x, y, σxxline, npts, R, θ0, L"\sigma_{xx}")
    PyPlot.subplot(3, 3, 8)
    circle_subplot(x, y, σyyline, npts, R, θ0, L"\sigma_{yy}")
    PyPlot.subplot(3, 3, 9)
    circle_subplot(x, y, σxyline, npts, R, θ0, L"\sigma_{xy}")
    PyPlot.suptitle(string(θ0), fontsize=30)
    PyPlot.tight_layout()

    # BEM solutions
    PyPlot.figure(figsize=(20, 20))
    PyPlot.subplot(3, 3, 1)
    circle_subplot(x, y, stresstrac[:, 1], npts, R, θ0, L"\sigma_{xx}")
    PyPlot.subplot(3, 3, 2)
    circle_subplot(x, y, stresstrac[:, 2], npts, R, θ0, L"\sigma_{yy}")
    PyPlot.subplot(3, 3, 3)
    circle_subplot(x, y, stresstrac[:, 3], npts, R, θ0, L"\sigma_{xy}")
    PyPlot.subplot(3, 3, 4)
    circle_subplot(x, y, stressdisp[:, 1], npts, R, θ0, L"\sigma_{xx}")
    PyPlot.subplot(3, 3, 5)
    circle_subplot(x, y, stressdisp[:, 2], npts, R, θ0, L"\sigma_{yy}")
    PyPlot.subplot(3, 3, 6)
    circle_subplot(x, y, stressdisp[:, 3], npts, R, θ0, L"\sigma_{xy}")
    PyPlot.subplot(3, 3, 7)
    circle_subplot(x, y, stresstrac[:, 1] + stressdisp[:, 1], npts, R, θ0, L"\sigma_{xx}")
    PyPlot.subplot(3, 3, 8)
    circle_subplot(x, y, stresstrac[:, 2] + stressdisp[:, 2], npts, R, θ0, L"\sigma_{yy}")
    PyPlot.subplot(3, 3, 9)
    circle_subplot(x, y, stresstrac[:, 3] + stressdisp[:, 3], npts, R, θ0, L"\sigma_{xy}")
    PyPlot.suptitle(string(θ0), fontsize=30)
    PyPlot.tight_layout()
    PyPlot.show()
end
ex_cylinder()
