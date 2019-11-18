using Revise
using PyCall
using PyPlot
using Bem2d

function plotlocal(els, idx, ufault, ufreesurfflat, ufreesurftopo, σfault, σfreesurfflat, σfreesurftopo, xobs, yobs, npts)
    # Pretty picture of displacements and stresses
    xfreesurf = unique([els.x1[idx["freesurfflat"]] ; els.x2[idx["freesurfflat"]]])
    xfill = [xfreesurf ; [10e3 ; -10e3 ; -10e3]]
    yfreesurf = unique([els.y1[idx["freesurfflat"]] ; els.y2[idx["freesurfflat"]]])
    yfill = [yfreesurf ; [5e3 ; 5e3 ; minimum(yfreesurf)]]

    # Combine two fields for total displacement and stress fields
    ufield = sqrt.((ufault - ufreesurfflat)[:, 1].^2 + (ufault - ufreesurfflat)[:, 2].^2)
    σxx = (σfault - σfreesurfflat)[:, 1]
    σyy = (σfault - σfreesurfflat)[:, 2]
    σxy = (σfault - σfreesurfflat)[:, 3]
    I1 = σxx + σyy  # 1st invariant
    I2 = σxx .* σyy - σxy.^2  # 2nd invariant
    J2 = (I1.^2) ./ 3.0 - I2  # 2nd invariant (deviatoric)
    σfield = log10.(abs.(J2))

    ncontours = 20
    ucontours = collect(LinRange(0.0, 1.0, ncontours))
    σcontours = collect(LinRange(10.0, 16.0, ncontours))

    figure(figsize = (15, 10))
    subplot(2, 2, 1)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ucontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$||u||$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ucontours, linewidths = 0.5, colors = "gray")
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
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, cmap = get_cmap("hot_r"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, linewidths = 0.5, colors = "gray")
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
    ufield = sqrt.((ufault - ufreesurftopo)[:, 1].^2 + (ufault - ufreesurftopo)[:, 2].^2)
    σxx = (σfault - σfreesurftopo)[:, 1]
    σyy = (σfault - σfreesurftopo)[:, 2]
    σxy = (σfault - σfreesurftopo)[:, 3]
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
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, cmap = get_cmap("hot_r"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), σcontours, linewidths = 0.5, colors = "gray")
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

function ex_thrusttopo()
    close("all")
    μ = 30e9
    ν = 0.25
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
    ∂const = init∂(els)
    ∂quad = init∂(els)

    #
    # CS elements
    #
    # Partial derivatves for solving BEM problem
    _, _, ∂const["t"]["fault"]["freesurfflat"] = ∂constuσ(slip2uσ, els, idx["fault"], idx["freesurfflat"], μ, ν)
    _, _, ∂const["t"]["freesurfflat"]["freesurfflat"] = ∂constuσ(slip2uσ, els, idx["freesurfflat"], idx["freesurfflat"], μ, ν)
    _, _, ∂const["t"]["fault"]["freesurftopo"] = ∂constuσ(slip2uσ, els, idx["fault"], idx["freesurftopo"], μ, ν)
    _, _, ∂const["t"]["freesurftopo"]["freesurftopo"] = ∂constuσ(slip2uσ, els, idx["freesurftopo"], idx["freesurftopo"], μ, ν)
    
    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(2 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    ufreesurfflatslip = inv(∂const["t"]["freesurfflat"]["freesurfflat"]) * ∂const["t"]["fault"]["freesurfflat"] * faultslip
    ufreesurftoposlip = inv(∂const["t"]["freesurftopo"]["freesurftopo"]) * ∂const["t"]["fault"]["freesurftopo"] * faultslip

    # Forward model for volumetric displacements and stresses
    ufault, σfault = constuσ(slip2uσ, xobs, yobs, els, idx["fault"], ones(size(idx["fault"])), zeros(size(idx["fault"])), μ, ν)
    ufreesurfflat, σfreesurfflat = constuσ(slip2uσ, xobs, yobs, els, idx["freesurfflat"], ufreesurfflatslip[1:2:end], ufreesurfflatslip[2:2:end], μ, ν)
    ufreesurftopo, σfreesurftopo = constuσ(slip2uσ, xobs, yobs, els, idx["freesurftopo"], ufreesurftoposlip[1:2:end], ufreesurftoposlip[2:2:end], μ, ν)
    
    # Pretty plotting
    plotlocal(els, idx, ufault, ufreesurfflat, ufreesurftopo, σfault, σfreesurfflat, σfreesurftopo, xobs, yobs, npts)

    #
    # 3QN elements
    #
    # Partial derivatves for solving BEM problem
    _ , _, ∂quad["t"]["fault"]["freesurfflat"] = ∂quaduσ(slip2uσ, els, idx["fault"], idx["freesurfflat"], μ, ν)
    _, _, ∂quad["t"]["freesurfflat"]["freesurfflat"] = ∂quaduσ(slip2uσ, els, idx["freesurfflat"], idx["freesurfflat"], μ, ν)
    _, _, ∂quad["t"]["fault"]["freesurftopo"] = ∂quaduσ(slip2uσ, els, idx["fault"], idx["freesurftopo"], μ, ν)
    _, _, ∂quad["t"]["freesurftopo"]["freesurftopo"] = ∂quaduσ(slip2uσ, els, idx["freesurftopo"], idx["freesurftopo"], μ, ν)

    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(6 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    ufreesurfflatslip = inv(∂quad["t"]["freesurfflat"]["freesurfflat"]) * ∂quad["t"]["fault"]["freesurfflat"] * faultslip
    ufreesurftoposlip = inv(∂quad["t"]["freesurftopo"]["freesurftopo"]) * ∂quad["t"]["fault"]["freesurftopo"] * faultslip

    # Forward model for volumetric displacements and stresses
    ufault, σfault = quaduσ(slip2uσ, xobs, yobs, els, idx["fault"], quadstack(ones(3 * length(idx["fault"]))), quadstack(zeros(3 * length(idx["fault"]))), μ, ν)
    ufreesurfflat, σfreesurfflat = quaduσ(slip2uσ, xobs, yobs, els, idx["freesurfflat"], quadstack(ufreesurfflatslip[1:2:end]), quadstack(ufreesurfflatslip[2:2:end]), μ, ν)
    ufreesurftopo, σfreesurftopo = quaduσ(slip2uσ, xobs, yobs, els, idx["freesurftopo"], quadstack(ufreesurftoposlip[1:2:end]), quadstack(ufreesurftoposlip[2:2:end]), μ, ν)

    # Pretty plotting
    plotlocal(els, idx, ufault, ufreesurfflat, ufreesurftopo, σfault, σfreesurfflat, σfreesurftopo, xobs, yobs, npts)
    
end
ex_thrusttopo()
