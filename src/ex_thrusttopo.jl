using Revise
using PyCall
using PyPlot
using Bem2d

function ex_thrusttopo()
    μ = 30e9
    ν = 0.25
    els = Elements()

    # Observation points for internal evaluation and visualization
    npts = 100
    xobs = LinRange(-10e3, 10e3, npts)
    yobs = LinRange(-5e3, 5e3, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # Topographic free surface
    nfreesurf = 20
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, nfreesurf)
    y1 = -1e3 * @.atan(x1 / 1e3)
    y2 = -1e3 * @.atan(x2 / 1e3)
    els.x1[els.endidx + 1 : els.endidx + nfreesurf] = x1
    els.y1[els.endidx + 1 : els.endidx + nfreesurf] = y1
    els.x2[els.endidx + 1 : els.endidx + nfreesurf] = x2
    els.y2[els.endidx + 1 : els.endidx + nfreesurf] = y2
    els.name[els.endidx + 1 : els.endidx + nfreesurf] .= "freesurf"
    standardize_elements!(els)

    # Curved fault
    nfault = 10
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, nfault)
    y1 = 3e3 * @.atan(x1 / 1e3)
    y2 = 3e3 * @.atan(x2 / 1e3)
    els.x1[els.endidx + 1 : els.endidx + nfault] = x1
    els.y1[els.endidx + 1 : els.endidx + nfault] = y1
    els.x2[els.endidx + 1 : els.endidx + nfault] = x2
    els.y2[els.endidx + 1 : els.endidx + nfault] = y2
    els.name[els.endidx + 1 : els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Partial derivatves
    srcidx = findall(x -> x == "fault", els.name)
    obsidx = findall(x -> x == "freesurf", els.name)
    d1, s1, t1 = ∂constslip(els, srcidx, obsidx, μ, ν)
    srcidx = findall(x -> x == "freesurf", els.name)
    obsidx = findall(x -> x == "freesurf", els.name)
    d2, s2, t2 = ∂constslip(els, srcidx, obsidx, μ, ν)

    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(2 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    ufreesurf = inv(t2) * t1 * faultslip

    # Fault in full space
    ufault = zeros(length(xobs), 2)
    σfault = zeros(length(xobs), 3)
    faultidx = findall(x -> x == "fault", els.name)
    for i in 1:length(faultidx)
        u, σ = constslip(xobs, yobs, els.halflength[faultidx[i]],
            μ, ν, 1, 0, els.xcenter[faultidx[i]], els.ycenter[faultidx[i]],
            els.rotmat[faultidx[i], :, :], els.rotmatinv[faultidx[i], :, :])
        ufault += u
        σfault += σ
    end

    # Free surface in full space
    ufreesurf = zeros(length(xobs), 2)
    σfreesurf = zeros(length(xobs), 3)
    freesurfidx = findall(x -> x == "freesurf", els.name)
    for i in 1:length(freesurfidx)
        u, σ = constslip(xobs, yobs, els.halflength[freesurfidx[i]],
            μ, ν, ufreesurf[1:2:end][i], ufreesurf[2:2:end][i],
            els.xcenter[freesurfidx[i]], els.ycenter[freesurfidx[i]],
            els.rotmat[freesurfidx[i], :, :], els.rotmatinv[freesurfidx[i], :, :])
        ufreesurf += u
        σfreesurf += σ
    end

    # Pretty of displacements and stresses
    freesurfidx = findall(x -> x == "freesurf", els.name)
    xfreesurf = unique([els.x1[freesurfidx] ; els.x2[freesurfidx]])
    xfill = [xfreesurf ; [10e3 ; -10e3 ; -10e3]]
    yfreesurf = unique([els.y1[freesurfidx] ; els.y2[freesurfidx]])
    yfill = [yfreesurf ; [5e3 ; 5e3 ; minimum(yfreesurf)]]
    ufield = @.sqrt((ufault + ufreesurf)[:, 1].^2 + (ufault + ufreesurf)[:, 2].^2)
    σxx = (σfreesurf + σfault)[:, 1]
    σyy = (σfreesurf + σfault)[:, 2]
    σxy = (σfreesurf + σfault)[:, 3]
    I1 = σxx + σyy  # 1st invariant
    I2 = σxx .* σyy - σxy.^2  # 2nd invariant
    J2 = (I1.^2) ./ 3.0 - I2  # 2nd invariant (deviatoric)
    σfield = @.log10(abs(J2))

    ncontours = 10
    figure(figsize=(6, 8))
    subplot(2, 1, 1)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ncontours, cmap=get_cmap("plasma"))
    colorbar(fraction=0.020, pad=0.05, extend="both", label=L"$||u||$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ncontours, linewidths=0.5, colors="gray")
    fill(xfill, yfill, "w", zorder=30)
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")
    title("displacement magnitude")

    subplot(2, 1, 2)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), ncontours, cmap=get_cmap("hot_r"))
    colorbar(fraction=0.020, pad=0.05, extend="both", label=L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), ncontours, linewidths=0.5, colors="gray")
    fill(xfill, yfill, "w", zorder=30)
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")
    title("2nd stress invariant (deviatoric)")
    show()
end
ex_thrusttopo()
