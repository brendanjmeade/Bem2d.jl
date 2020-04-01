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

function circle_subplot(nrows, ncols, plotidx, x, y, mat, npts, R, θ0, title_string)
    fontsize = 20
    contour_levels = 20
    contour_color = "white"
    contour_linewidth = 0.5
    color_scale = 1

    subplot(nrows, ncols, plotidx)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels)
    cbar = colorbar(fraction=0.020, pad=0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    contour(reshape(x, npts, npts), reshape(y, npts, npts), reshape(mat, npts, npts), levels=contour_levels, colors=contour_color, linewidths=contour_linewidth)
    title(title_string, fontsize=fontsize)

    # Draw entire circle and region of applied tractions
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, 360)
    plot([x1, x2], [y1, y2], "-k", linewidth=2)
    x1, y1, x2, y2 = discretized_arc(-θ0, θ0, R, 50)
    plot([x1, x2], [y1, y2], "-r", linewidth=2)
    x1, y1, x2, y2 = discretized_arc(-θ0+deg2rad(180), θ0+deg2rad(180), R, 50)
    plot([x1, x2], [y1, y2], "-r", linewidth=2)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)
end

function ex_braziltest()
    PyPlot.close("all")
    mu = 3e10
    nu = 0.25

    # Solution from Hondros (1959) as summarized by Wei and Chau 2013
    p = -1.0e5 # Applied radial pressure over arc
    theta0 = deg2rad(1.0) # Arc length over which pressure is applied
    nels = 360
    R = nels / (2 * pi) # Radius of disc ensuring that all elements are 1 unit long
    mmax = 1000 # Max number of terms in Hondros series
    npts = 50
    x, y = Bem2d.obsgrid(-R, -R, R, R, npts)
    r = @. sqrt(x^2 + y^2)
    θ = @. atan(y, x)

    #! Analytic stresses in cylindrical coordinates
    Srr = zeros(length(x))
    Sθθ = zeros(length(x))
    Srθ = zeros(length(x))
    Srrconstterm = 2.0 * theta0 * -p / deg2rad(180)
    Sθθconstterm = 2.0 * theta0 * -p / deg2rad(180)
    leadingterm = 2.0 * -p / deg2rad(180)
    for m in 1:mmax
        Srr += @. (r/R)^(2*m-2) * (1-(1-1/m)*(r/R)^2) * sin(2*m*theta0) * cos(2*m*θ)
        Sθθ += @. (r/R)^(2*m-2) * (1-(1+1/m)*(r/R)^2) * sin(2*m*theta0) * cos(2*m*θ)
        Srθ += @. ((r/R)^(2*m) - (r/R)^(2*m-2)) * sin(2*m*theta0) * sin(2*m*θ)
    end
    Srr = @. Srrconstterm + leadingterm * Srr
    Sθθ = @. Sθθconstterm - leadingterm * Sθθ
    Srθ = @. leadingterm * Srθ

    #! Convert analytic cylindrical stresses to Cartesian
    Sxx = zeros(length(x))
    Syy = zeros(length(x))
    Sxy = zeros(length(x))
    for i in 1:length(x) # Project a single matrix to Cartesian coordinates
        cylindrical_stress_tensor = [Srr[i] Srθ[i] ; Srθ[i] Sθθ[i]]
        transformation_matrix = [cos(θ[i]) -sin(θ[i]) ; sin(θ[i]) cos(θ[i])]
        cartesian_stress_tensor = transformation_matrix * cylindrical_stress_tensor * transpose(transformation_matrix)
        Sxx[i] = cartesian_stress_tensor[1, 1]
        Syy[i] = cartesian_stress_tensor[2, 2]
        Sxy[i] = cartesian_stress_tensor[1, 2]
    end

    #! BEM solution
    els = Bem2d.Elements(Int(1e5))
    x1, y1, x2, y2 = discretized_arc(deg2rad(-180), deg2rad(180), R, nels)
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
        ytrac[i] = -ytrac[i] # TODO.  FIX ME.  Just a desperate attempt
    end

    #! Zero out the tractions on the area without contact
    deleteidx = findall(x -> (x>theta0 && x<deg2rad(180)-theta0), thetaels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0
    deleteidx = findall(x -> (x<-theta0 && x>-deg2rad(180)+theta0), thetaels)
    xtrac[deleteidx] .= 0
    ytrac[deleteidx] .= 0

    #! Scale tractions by element lengths
    xtracscaled = xtrac ./ els.length[1:1:els.endidx]
    ytracscaled = ytrac ./ els.length[1:1:els.endidx]

    #! Kernels
    # U*: traction to displacement
    # T*: displacement to displacement
    # A*: traction to traction
    # H*: displacement to traction
    Tstar, Sstar, Hstar = partialsconstdispstress(slip2dispstress, els, idx["circle"], idx["circle"], mu, nu)
    Ustar, Dstar, Astar = partialsconstdispstress(trac2dispstress, els, idx["circle"], idx["circle"], mu, nu)

    #! Induced displacements from Direct BEM
    dispall = (inv(Tstar + 0.5 * I(size(Tstar)[1]))) * Ustar * interleave(xtrac, ytrac)

    #! Displacement discontinuity method (indirect)
    dispall = Ustar * interleave(xtrac, ytrac)

    #! Streses from tractions
    Dtrac, Strac = constdispstress(trac2dispstress, x, y, els, idx["circle"], xtracscaled, ytracscaled, mu, nu)

    #! Stresses from traction induced displacements
    Ddisp, Sdisp = constdispstress(slip2dispstress, x, y, els, idx["circle"], dispall[1:2:end], dispall[2:2:end], mu, nu)

    # #! Isolate the values inside the circle
    to_nan_idx = findall(x -> x > 0.9 * R, r)
    Srr[to_nan_idx] .= NaN
    Sθθ[to_nan_idx] .= NaN
    Srθ[to_nan_idx] .= NaN
    Sxx[to_nan_idx] .= NaN
    Syy[to_nan_idx] .= NaN
    Sxy[to_nan_idx] .= NaN
    Strac[to_nan_idx, 1] .= NaN
    Strac[to_nan_idx, 2] .= NaN
    Strac[to_nan_idx, 3] .= NaN
    Sdisp[to_nan_idx, 1] .= NaN
    Sdisp[to_nan_idx, 2] .= NaN
    Sdisp[to_nan_idx, 3] .= NaN
    Sbem = @. Strac + Sdisp

    #! Summary figure
    fontsize = 24
    figure(figsize=(30,20))
    nrows = 3
    ncols = 6

    #! Traction forcing
    subplot(3, 6, 1)
    plot([els.x1[1:els.endidx], els.x2[1:els.endidx]], [els.y1[1:els.endidx], els.y2[1:els.endidx]], "-k") 
    quiver(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], xtrac, ytrac)
    title("applied tractions", fontsize=fontsize)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    #! Induced displacements
    subplot(3, 6, 2)
    plot([els.x1[1:els.endidx], els.x2[1:els.endidx]], [els.y1[1:els.endidx], els.y2[1:els.endidx]], "-k") 
    quiver(els.xcenter[1:1:els.endidx], els.ycenter[1:1:els.endidx], dispall[1:2:end], dispall[2:2:end], color="green")
    title("induced displacements", fontsize=fontsize)
    gca().set_aspect("equal")
    gca().tick_params(labelsize=fontsize)

    #! Analytic solutions
    circle_subplot(nrows, ncols, 7, x, y, Srr, npts, R, theta0, "Srr (analytic}")
    circle_subplot(nrows, ncols, 8, x, y, Sθθ, npts, R, theta0, "Sthetatheta (analytic)")
    circle_subplot(nrows, ncols, 9, x, y, Srθ, npts, R, theta0, "Srtheta(analytic)")
    circle_subplot(nrows, ncols, 13, x, y, Sxx, npts, R, theta0, "Sxx (analytic)")
    circle_subplot(nrows, ncols, 14, x, y, Syy, npts, R, theta0, "Syy (analytic)")
    circle_subplot(nrows, ncols, 15, x, y, Sxy, npts, R, theta0, "Sxy (analytic)")

    #! BEM solutions
    circle_subplot(nrows, ncols, 4, x, y, Strac[:, 1], npts, R, theta0, "Sxx (tractions only)")
    circle_subplot(nrows, ncols, 5, x, y, Strac[:, 2], npts, R, theta0, "Syy (tractions only)")
    circle_subplot(nrows, ncols, 6, x, y, Strac[:, 3], npts, R, theta0, "Sxy (tractions only)")
    circle_subplot(nrows, ncols, 10, x, y, Sdisp[:, 1], npts, R, theta0, "Sxx (DDM)")
    circle_subplot(nrows, ncols, 11, x, y, Sdisp[:, 2], npts, R, theta0, "Syy (DDM)")
    circle_subplot(nrows, ncols, 12, x, y, Sdisp[:, 3], npts, R, theta0, "Sxy (DDM)")
    circle_subplot(nrows, ncols, 16, x, y, Sbem[:, 1], npts, R, theta0, "Sxx( sum)")
    circle_subplot(nrows, ncols, 17, x, y, Sbem[:, 2], npts, R, theta0, "Syy (sum)")
    circle_subplot(nrows, ncols, 18, x, y, Sbem[:, 3], npts, R, theta0, "Sxy (sum)")
    tight_layout()
    show()
end
ex_braziltest()
