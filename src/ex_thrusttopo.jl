using Revise
using PyCall
using PyPlot
using Bem2d

function ex_thrusttopo()
    mu = 30e9
    nu = 0.25
    els = Elements()

    # Observation points for internal evaluation and visualization
    npts = 100
    xobs = LinRange(-10e3, 10e3, npts)
    yobs = LinRange(-5e3, 5e3, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # Topographic free surface
    nfreesurface = 20
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, nfreesurface)
    y1 = -1e3 * @.atan(x1 / 1e3)
    y2 = -1e3 * @.atan(x2 / 1e3)
    els.x1[els.endidx + 1 : els.endidx + nfreesurface] = x1
    els.y1[els.endidx + 1 : els.endidx + nfreesurface] = y1
    els.x2[els.endidx + 1 : els.endidx + nfreesurface] = x2
    els.y2[els.endidx + 1 : els.endidx + nfreesurface] = y2
    els.name[els.endidx + 1 : els.endidx + nfreesurface] .= "freesurface"
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
    obsidx = findall(x -> x == "freesurface", els.name)
    d1, s1, t1 = partials_constslip(els, srcidx, obsidx, mu, nu)
    srcidx = findall(x -> x == "freesurface", els.name)
    obsidx = findall(x -> x == "freesurface", els.name)
    d2, s2, t2 = partials_constslip(els, srcidx, obsidx, mu, nu)

    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(2 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    ufreesurface = inv(t2) * t1 * faultslip

    # Fault in full space
    ufault = zeros(length(xobs), 2)
    σfault = zeros(length(xobs), 3)
    faultidx = findall(x -> x == "fault", els.name)
    for i in 1:length(faultidx)
        u, σ = dispstress_constslip(xobs, yobs, els.halflength[faultidx[i]],
            mu, nu, 1, 0, els.xcenter[faultidx[i]], els.ycenter[faultidx[i]],
            els.rotmat[faultidx[i], :, :], els.rotmatinv[faultidx[i], :, :])
        ufault += u
        σfault += σ
    end

    # Free surface in full space
    ufreesurface = zeros(length(xobs), 2)
    σfreesurface = zeros(length(xobs), 3)
    freesurfaceidx = findall(x -> x == "freesurface", els.name)
    for i in 1:length(freesurfaceidx)
        u, σ = dispstress_constslip(xobs, yobs, els.halflength[freesurfaceidx[i]],
            mu, nu, ufreesurface[1:2:end][i], ufreesurface[2:2:end][i],
            els.xcenter[freesurfaceidx[i]], els.ycenter[freesurfaceidx[i]],
            els.rotmat[freesurfaceidx[i], :, :], els.rotmatinv[freesurfaceidx[i], :, :])
        ufreesurface += u
        σfreesurface += σ
    end

    # Pretty of displacements and stresses
    freesurfaceidx = findall(x -> x == "freesurface", els.name)
    xfreesurface = unique([els.x1[freesurfaceidx] ; els.x2[freesurfaceidx]])
    xfill = [xfreesurface ; [10e3 ; -10e3 ; -10e3]]
    yfreesurface = unique([els.y1[freesurfaceidx] ; els.y2[freesurfaceidx]])
    yfill = [yfreesurface ; [5e3 ; 5e3 ; minimum(yfreesurface)]]
    ufield = @.sqrt((ufault + ufreesurface)[:, 1].^2 + (ufault + ufreesurface)[:, 2].^2)
    σxx = (σfreesurface + σfault)[:, 1]
    σyy = (σfreesurface + σfault)[:, 2]
    σxy = (σfreesurface + σfault)[:, 3]
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
