using Revise
using PyCall
using PyPlot
using Bem2d

function ex_thrusttopo()
    mu = 30e9
    nu = 0.25
    elements = Elements()

    # Observation points for internal evaluation and visualization
    npts= 200
    xobs = LinRange(-10e3, 10e3, npts)
    yobs = LinRange(-5e3, 5e3, npts)
    xobs, yobs = meshgrid(xobs, yobs)
    xobs = xobs[:]
    yobs = yobs[:]

    # Topographic free surface
    x1, y1, x2, y2 = discretizedline(-10e3, 0, 10e3, 0, 200)
    y1 = -1e3 * @.atan(x1 / 1e3)
    y2 = -1e3 * @.atan(x2 / 1e3)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "freesurface"
    end
    standardize_elements!(elements)

    # Curved fault
    x1, y1, x2, y2 = discretizedline(-7e3, 0e3, 0, 0, 100)
    y1 = 3e3 * @.atan(x1 / 1e3)
    y2 = 3e3 * @.atan(x2 / 1e3)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "fault"
    end
    standardize_elements!(elements)

    # Partial derivatves
    srcidx = findall(x -> x == "fault", elements.name)
    obsidx = findall(x -> x == "freesurface", elements.name)
    d1, s1, t1 = partials_constslip(elements, srcidx, obsidx, mu, nu)
    srcidx = findall(x -> x == "freesurface", elements.name)
    obsidx = findall(x -> x == "freesurface", elements.name)
    d2, s2, t2 = partials_constslip(elements, srcidx, obsidx, mu, nu)

    # Remove and separate BCs, local to global transform here?
    nfaultelements = length(findall(x -> x == "fault", elements.name))
    faultslip = zeros(2 * nfaultelements)
    faultslip[1:2:end] .= 1.0 # Global coordinate system

    # Solve the BEM problem
    disp_freesurface = inv(t2) * t1 * faultslip

    # Fault in full space
    dispfault = zeros(length(xobs), 2)
    stressfault = zeros(length(xobs), 3)
    faultidx = findall(x -> x == "fault", elements.name)
    for i in 1:length(faultidx)
        disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[faultidx[i]],
            mu, nu, 1, 0, elements.xcenter[faultidx[i]], elements.ycenter[faultidx[i]],
            elements.rotmat[faultidx[i], :, :], elements.rotmatinv[faultidx[i], :, :])
        dispfault += disp
        stressfault += stress
    end

    # Free surface in full space
    dispfreesurface::Array{Float64} = zeros(length(xobs), 2)
    stressfreesurface::Array{Float64} = zeros(length(xobs), 3)
    freesurfaceidx::Array{Int64} = findall(x -> x == "freesurface", elements.name)
    for i in 1:length(freesurfaceidx)
        disp, stress = dispstress_constslip(xobs, yobs, elements.halflength[freesurfaceidx[i]],
            mu, nu, disp_freesurface[1:2:end][i], disp_freesurface[2:2:end][i],
            elements.xcenter[freesurfaceidx[i]], elements.ycenter[freesurfaceidx[i]],
            elements.rotmat[freesurfaceidx[i], :, :], elements.rotmatinv[freesurfaceidx[i], :, :])
        dispfreesurface += disp
        stressfreesurface += stress
    end

    # Pretty of displacements and stresses
    freesurfaceidx = findall(x -> x == "freesurface", elements.name)
    xfreesurface = unique([elements.x1[freesurfaceidx] ; elements.x2[freesurfaceidx]])
    xfill = [xfreesurface ; [10e3 ; -10e3 ; -10e3]]
    yfreesurface = unique([elements.y1[freesurfaceidx] ; elements.y2[freesurfaceidx]])
    yfill = [yfreesurface ; [5e3 ; 5e3 ; minimum(yfreesurface)]]
    ufield = @.sqrt((dispfault + dispfreesurface)[:, 1].^2 + (dispfault + dispfreesurface)[:, 2].^2)
    σxx = (stressfreesurface + stressfault)[:, 1]
    σyy = (stressfreesurface + stressfault)[:, 2]
    σxy = (stressfreesurface + stressfault)[:, 3]
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
    plotelements(elements)
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
    plotelements(elements)
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")
    title("2nd stress invariant (deviatoric)")
    show()
end
ex_thrusttopo()
