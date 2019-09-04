using Revise
using PyCall
using PyPlot
using Bem2d

function ex_thrusttopo()
    μ = 30e9
    ν = 0.25
    els = Elements()

    # Observation points for internal evaluation and visualization
    npts = 300
    xobs, yobs = obsgrid(-10e3, -5e3, 10e3, 5e3, npts)

    # Topographic free surface
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, nfreesurf)
    y1 = -1e3 * @.atan(x1 / 1e3)
    y2 = -1e3 * @.atan(x2 / 1e3)
    els.x1[els.endidx + 1:els.endidx + nfreesurf] = x1
    els.y1[els.endidx + 1:els.endidx + nfreesurf] = y1
    els.x2[els.endidx + 1:els.endidx + nfreesurf] = x2
    els.y2[els.endidx + 1:els.endidx + nfreesurf] = y2
    els.name[els.endidx + 1:els.endidx + nfreesurf] .= "freesurf"
    standardize_elements!(els)

    # Curved fault
    nfault = 100
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, nfault)
    y1 = 3e3 * @.atan(x1 / 1e3)
    y2 = 3e3 * @.atan(x2 / 1e3)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Partial derivatves
    srcidx = findall(x->x == "fault", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u1, ∂σ1, ∂t1 = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)
    srcidx = findall(x->x == "freesurf", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u2, ∂σ2, ∂t2 = ∂constuσ(slip2uσ, els, srcidx, obsidx, μ, ν)

    # Solve the BEM problem for unit slip in the x-direction
    faultslip = zeros(2 * nfault)
    faultslip[1:2:end] .= 1.0 # Global coordinate system
    ufreesurf = inv(∂t2) * ∂t1 * faultslip

    # Fault in full space
    faultidx = findall(x->x == "fault", els.name)
    ufault, σfault = constuσ(slip2uσ, xobs, yobs, els, faultidx,
        ones(size(faultidx)), zeros(size(faultidx)), μ, ν)    

    # Free surface in full space
    freesurfidx = findall(x->x == "faultsurf", els.name)
    ufreesurf, σfreesurf = constuσ(slip2uσ, xobs, yobs, els, freesurfidx,
        ufreesurf[1:2:end], ufreesurf[2:2:end], μ, ν)

    # Pretty picture of displacements and stresses
    freesurfidx = findall(x->x == "freesurf", els.name)
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
    figure(figsize = (8, 8))
    subplot(2, 1, 1)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ncontours, cmap = get_cmap("plasma"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$||u||$ (m)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(ufield, npts, npts), ncontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")

    subplot(2, 1, 2)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), ncontours, cmap = get_cmap("hot_r"))
    colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)")
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts),
        reshape(σfield, npts, npts), ncontours, linewidths = 0.5, colors = "gray")
    fill(xfill, yfill, "w", zorder = 30)
    plotelements(els)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")
    show()
end
ex_thrusttopo()
