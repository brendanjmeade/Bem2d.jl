using Revise
using PyCall
using PyPlot
using Bem2d

function asdf(els)
    idx = Dict()
    names = unique(els.name)
    for i in 1:length(names)
        idx[names[i]] = getidx(names[i], els)
    end
    return idx
end


function ex_thrusttopo()
    μ = 30e9
    ν = 0.25
    els = Elements(Int(1e5))

    # Observation points for internal evaluation and visualization
    npts = 50
    xobs, yobs = obsgrid(-10e3, -5e3, 10e3, 5e3, npts)
    nfreesurf = 40
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

    # Dictionary with indices to all names/labels
    idx = asdf(els)
    
    # Partial derivatves
    ∂u1, ∂σ1, ∂t1 = ∂constuσ(slip2uσ, els, idx["fault"], idx["freesurfflat"], μ, ν)
    ∂u2, ∂σ2, ∂t2 = ∂constuσ(slip2uσ, els, idx["freesurfflat"], idx["freesurfflat"], μ, ν)
    
    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(2 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    ufreesurfslip = inv(∂t2) * ∂t1 * faultslip

    #
    # Fault in half space
    #
    ufault, σfault = constuσ(slip2uσ, xobs, yobs, els, idx["fault"], ones(size(idx["fault"])), zeros(size(idx["fault"])), μ, ν)
    ufreesurf, σfreesurf = constuσ(slip2uσ, xobs, yobs, els, idx["freesurfflat"], ufreesurfslip[1:2:end], ufreesurfslip[2:2:end], μ, ν)

    # Pretty picture of displacements and stresses
    xfreesurf = unique([els.x1[idx["freesurfflat"]] ; els.x2[idx["freesurfflat"]]])
    xfill = [xfreesurf ; [10e3 ; -10e3 ; -10e3]]
    yfreesurf = unique([els.y1[idx["freesurfflat"]] ; els.y2[idx["freesurfflat"]]])
    yfill = [yfreesurf ; [5e3 ; 5e3 ; minimum(yfreesurf)]]

    # Combine two fields for total displacement and stress fields
    ufield = sqrt.((ufault - ufreesurf)[:, 1].^2 + (ufault - ufreesurf)[:, 2].^2)
    σxx = (σfault - σfreesurf)[:, 1]
    σyy = (σfault - σfreesurf)[:, 2]
    σxy = (σfault - σfreesurf)[:, 3]
    I1 = σxx + σyy  # 1st invariant
    I2 = σxx .* σyy - σxy.^2  # 2nd invariant
    J2 = (I1.^2) ./ 3.0 - I2  # 2nd invariant (deviatoric)
    σfield = log10.(abs.(J2))

    ncontours = 20
    ucontours = collect(LinRange(0.0, 1.0, ncontours))
    σcontours = collect(LinRange(10.0, 16.0, ncontours))

    close("all")
    figure(figsize = (15, 10))
    subplot(2, 2, 1)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ucontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$||u||$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ucontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 2, 3)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, cmap = get_cmap("hot_r"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")


    #
    # Topography
    #
    ∂u1, ∂σ1, ∂t1 = ∂constuσ(slip2uσ, els, idx["fault"], idx["freesurftopo"], μ, ν)
    ∂u2, ∂σ2, ∂t2 = ∂constuσ(slip2uσ, els, idx["freesurftopo"], idx["freesurftopo"], μ, ν)

    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(2 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    ufreesurfslip = inv(∂t2) * ∂t1 * faultslip

    # Fault in full space
    ufault, σfault = constuσ(slip2uσ, xobs, yobs, els, idx["fault"], ones(size(idx["fault"])), zeros(size(idx["fault"])), μ, ν)

    # Free surface in full space
    ufreesurf, σfreesurf = constuσ(slip2uσ, xobs, yobs, els, idx["freesurftopo"], ufreesurfslip[1:2:end], ufreesurfslip[2:2:end], μ, ν)

    # Pretty picture of displacements and stresses
    xfreesurf = unique([els.x1[idx["freesurftopo"]] ; els.x2[idx["freesurftopo"]]])
    xfill = [xfreesurf ; [10e3 ; -10e3 ; -10e3]]
    yfreesurf = unique([els.y1[idx["freesurftopo"]] ; els.y2[idx["freesurftopo"]]])
    yfill = [yfreesurf ; [5e3 ; 5e3 ; minimum(yfreesurf)]]

    # Combine two fields for total displacement and stress fields
    ufield = sqrt.((ufault - ufreesurf)[:, 1].^2 + (ufault - ufreesurf)[:, 2].^2)
    σxx = (σfault - σfreesurf)[:, 1]
    σyy = (σfault - σfreesurf)[:, 2]
    σxy = (σfault - σfreesurf)[:, 3]
    I1 = σxx + σyy  # 1st invariant
    I2 = σxx .* σyy - σxy.^2  # 2nd invariant
    J2 = (I1.^2) ./ 3.0 - I2  # 2nd invariant (deviatoric)
    σfield = log10.(abs.(J2))

    subplot(2, 2, 2)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ucontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$||u||$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ucontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 2, 4)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, cmap = get_cmap("hot_r"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    for i in 1:els.endidx
        plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40)
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
end
ex_thrusttopo()
