using Revise
using PyCall
using PyPlot
using Bem2d

"""
    plotlocal(els, idx, dispfault, dispfreesurfflat, dispfreesurftopo, stressfault, stressfreesurfflat, stressfreesurftopo, xobs, yobs, npts)

Generate four panel plots to shw displacements and stresses for the
flat and topographic cases.
"""
function plotlocal(els, idx, dispfault, dispfreesurfflat, dispfreesurftopo, stressfault, stressfreesurfflat, stressfreesurftopo, xobs, yobs, npts)
    # Pretty picture of displacements and stresses
    xfreesurf = unique([els.x1[idx["freesurfflat"]] ; els.x2[idx["freesurfflat"]]])
    xfill = [xfreesurf ; [10e3 ; -10e3 ; -10e3]]
    yfreesurf = unique([els.y1[idx["freesurfflat"]] ; els.y2[idx["freesurfflat"]]])
    yfill = [yfreesurf ; [5e3 ; 5e3 ; minimum(yfreesurf)]]

    # Combine two fields for total displacement and stress fields
    dispfield = sqrt.((dispfault - dispfreesurfflat)[:, 1].^2 + (dispfault - dispfreesurfflat)[:, 2].^2)
    stressxx = (stressfault - stressfreesurfflat)[:, 1]
    stressyy = (stressfault - stressfreesurfflat)[:, 2]
    stressxy = (stressfault - stressfreesurfflat)[:, 3]
    I1 = stressxx + stressyy  # 1st invariant
    I2 = stressxx .* stressyy - stressxy.^2  # 2nd invariant
    J2 = (I1.^2) ./ 3.0 - I2  # 2nd invariant (deviatoric)
    stressfield = log10.(abs.(J2))

    ncontours = 20
    dispcontours = collect(LinRange(0.0, 1.0, ncontours))
    stresscontours = collect(LinRange(10.0, 16.0, ncontours))

    figure(figsize = (15, 10))
    subplot(2, 2, 1)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(dispfield, npts, npts), dispcontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$||u||$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(dispfield, npts, npts), dispcontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        if els.name[i] == "fault" || els.name[i] == "freesurfflat"
            plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
        end
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 2, 3)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(stressfield, npts, npts), stresscontours, cmap = get_cmap("hot_r"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(stressfield, npts, npts), stresscontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        if els.name[i] == "fault" || els.name[i] == "freesurfflat"
            plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
        end
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    # Topography
    xfreesurf = unique([els.x1[idx["freesurftopo"]] ; els.x2[idx["freesurftopo"]]])
    xfill = [xfreesurf ; [10e3 ; -10e3 ; -10e3]]
    yfreesurf = unique([els.y1[idx["freesurftopo"]] ; els.y2[idx["freesurftopo"]]])
    yfill = [yfreesurf ; [5e3 ; 5e3 ; minimum(yfreesurf)]]

    # Combine two fields for total displacement and stress fields
    dispfield = sqrt.((dispfault - dispfreesurftopo)[:, 1].^2 + (dispfault - dispfreesurftopo)[:, 2].^2)
    stressxx = (stressfault - stressfreesurftopo)[:, 1]
    stressyy = (stressfault - stressfreesurftopo)[:, 2]
    stressxy = (stressfault - stressfreesurftopo)[:, 3]
    I1 = stressxx + stressyy  # 1st invariant
    I2 = stressxx .* stressyy - stressxy.^2  # 2nd invariant
    J2 = (I1.^2) ./ 3.0 - I2  # 2nd invariant (deviatoric)
    stressfield = log10.(abs.(J2))

    subplot(2, 2, 2)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(dispfield, npts, npts), dispcontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$||u||$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(dispfield, npts, npts), dispcontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        if els.name[i] == "fault" || els.name[i] == "freesurftopo"
            plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
        end
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 2, 4)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(stressfield, npts, npts), stresscontours, cmap = get_cmap("hot_r"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(stressfield, npts, npts), stresscontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        if els.name[i] == "fault" || els.name[i] == "freesurftopo"
            plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
        end
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")
    tight_layout()
    show()
    return nothing
end

"""
    thrustfaulttopography()

Comparison of volume displacemnts and stresses from a fault with
and without topography.  Includes both constant and quadratic
element cases.
"""
function thrustfaulttopography()
    close("all")
    mu = 30e9
    nu = 0.25
    els = Elements(Int(1e5))

    # Observation points for internal evaluation and visualization
    npts = 50
    xobs, yobs = obsgrid(-10e3, -5e3, 10e3, 5e3, npts)
    nfreesurf = 20
    nfault = 10
    
    # Flat free surface
    x1, y1, x2, y2 = discretizedline(-20e3, 0, 20e3, 0, nfreesurf)
    y1 = -1e-3 * atan.(x1 / 1e3)
    y2 = -1e-3 * atan.(x2 / 1e3)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "freesurfflat"
    end
    standardize_elements!(els)

    # Topographic free surface
    x1, y1, x2, y2 = discretizedline(-20e3, 0, 20e3, 0, nfreesurf)
    y1 = -1e3 * atan.(x1 / 1e3)
    y2 = -1e3 * atan.(x2 / 1e3)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "freesurftopo"
    end
    standardize_elements!(els)
    
    # Curved fault
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, nfault)
    y1 = 3e3 * atan.(x1 / 1e3)
    y2 = 3e3 * atan.(x2 / 1e3)
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    standardize_elements!(els)

    # Create convience tools
    idx = getidxdict(els)
    partialsconst = initpartials(els)
    partialsquad = initpartials(els)

    #! CS elements
    # Partial derivatves for solving BEM problem
    _, _, partialsconst["trac"]["fault"]["freesurfflat"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurfflat"], mu, nu)
    _, _, partialsconst["trac"]["freesurfflat"]["freesurfflat"] = partialsconstdispstress(slip2dispstress, els, idx["freesurfflat"], idx["freesurfflat"], mu, nu)
    _, _, partialsconst["trac"]["fault"]["freesurftopo"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["freesurftopo"], mu, nu)
    _, _, partialsconst["trac"]["freesurftopo"]["freesurftopo"] = partialsconstdispstress(slip2dispstress, els, idx["freesurftopo"], idx["freesurftopo"], mu, nu)
    
    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(2 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    dispfreesurfflatslip = inv(partialsconst["trac"]["freesurfflat"]["freesurfflat"]) * partialsconst["trac"]["fault"]["freesurfflat"] * faultslip
    dispfreesurftoposlip = inv(partialsconst["trac"]["freesurftopo"]["freesurftopo"]) * partialsconst["trac"]["fault"]["freesurftopo"] * faultslip

    # Forward model for volumetric displacements and stresses
    dispfault, stressfault = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], ones(size(idx["fault"])), zeros(size(idx["fault"])), mu, nu)
    dispfreesurfflat, stressfreesurfflat = constdispstress(slip2dispstress, xobs, yobs, els, idx["freesurfflat"], dispfreesurfflatslip[1:2:end], dispfreesurfflatslip[2:2:end], mu, nu)
    dispfreesurftopo, stressfreesurftopo = constdispstress(slip2dispstress, xobs, yobs, els, idx["freesurftopo"], dispfreesurftoposlip[1:2:end], dispfreesurftoposlip[2:2:end], mu, nu)
    
    # Pretty plotting
    plotlocal(els, idx, dispfault, dispfreesurfflat, dispfreesurftopo, stressfault, stressfreesurfflat, stressfreesurftopo, xobs, yobs, npts)

    #! 3QN elements
    # Partial derivatves for solving BEM problem
    _, _, partialsquad["trac"]["fault"]["freesurfflat"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["freesurfflat"], mu, nu)
    _, _, partialsquad["trac"]["freesurfflat"]["freesurfflat"] = partialsquaddispstress(slip2dispstress, els, idx["freesurfflat"], idx["freesurfflat"], mu, nu)
    _, _, partialsquad["trac"]["fault"]["freesurftopo"] = partialsquaddispstress(slip2dispstress, els, idx["fault"], idx["freesurftopo"], mu, nu)
    _, _, partialsquad["trac"]["freesurftopo"]["freesurftopo"] = partialsquaddispstress(slip2dispstress, els, idx["freesurftopo"], idx["freesurftopo"], mu, nu)

    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(6 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    dispfreesurfflatslip = inv(partialsquad["trac"]["freesurfflat"]["freesurfflat"]) * partialsquad["trac"]["fault"]["freesurfflat"] * faultslip
    dispfreesurftoposlip = inv(partialsquad["trac"]["freesurftopo"]["freesurftopo"]) * partialsquad["trac"]["fault"]["freesurftopo"] * faultslip

    # Forward model for volumetric displacements and stresses
    dispfault, stressfault = quaddispstress(slip2dispstress, xobs, yobs, els, idx["fault"], quadstack(ones(3 * length(idx["fault"]))), quadstack(zeros(3 * length(idx["fault"]))), mu, nu)
    dispfreesurfflat, stressfreesurfflat = quaddispstress(slip2dispstress, xobs, yobs, els, idx["freesurfflat"], quadstack(dispfreesurfflatslip[1:2:end]), quadstack(dispfreesurfflatslip[2:2:end]), mu, nu)
    dispfreesurftopo, stressfreesurftopo = quaddispstress(slip2dispstress, xobs, yobs, els, idx["freesurftopo"], quadstack(dispfreesurftoposlip[1:2:end]), quadstack(dispfreesurftoposlip[2:2:end]), mu, nu)

    # Pretty plotting
    plotlocal(els, idx, dispfault, dispfreesurfflat, dispfreesurftopo, stressfault, stressfreesurfflat, stressfreesurftopo, xobs, yobs, npts)
    
end
thrustfaulttopography()
