using PyPlot
using JLD2
using FileIO
using DifferentialEquations
using Printf
using Infiltrator
using Bem2d

basename = "2020-02-10T13:32:04.689"

# Convert an ordered list of .png frames to a .mp4 animation
function png2mp4(framerate, outfilename)
    outfilename = outfilename * ".mp4"
    ffmpegstring = "ffmpeg -pattern_type glob -i '*.png' -c:v libx264 -r " * framerate * " -pix_fmt yuv420p " * outfilename
    run(`$ffmpegstring`)
end

function plotlocal(els, idx, dispfault, stressfault, xobs, yobs, npts, figuretitle)
    fontsize = 24
    ncontours = 20
    dispcontours = collect(-14:0.1:1)
    dispcontours = collect(LinRange(0, 10, 9))
    stresscontours = collect(-18:2:18)

    # Combine two fields for total displacement and stress fields
    dispfield = sqrt.(dispfault[:, 1].^2 + dispfault[:, 2].^2)
    stressxx = stressfault[:, 1]
    stressyy = stressfault[:, 2]
    stressxy = stressfault[:, 3]
    I1 = stressxx + stressyy  # 1st invariant
    I2 = stressxx .* stressyy - stressxy.^2  # 2nd invariant
    J2 = (I1.^2) ./ 3.0 - I2  # 2nd invariant (deviatoric)
    stressfield = log10.(abs.(J2))

    # Start figure generation
    figure(figsize = (12, 10))
    ax = subplot(2, 1, 1)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(dispfield, npts, npts), dispcontours, cmap = get_cmap("cool"))
    cbar = colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.set_label(label = L"$||v||$ (m)", fontsize=fontsize)
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(dispfield, npts, npts), dispcontours, linewidths = 0.5, colors = "black")
    for i in 1:els.endidx
        # if els.name[i] == "fault" | "freesurftopo"
            plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40, color="lightgray")
        # end
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    ylabel(L"$y$ (m)", fontsize=fontsize)
    ax.tick_params(labelsize = fontsize)

    ax = subplot(2, 1, 2)
    contourf(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(stressfield, npts, npts), stresscontours, cmap = rycroftcmap())
    cbar = colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.set_label(label = L"$\log_{10}|\mathrm{J}_2|$ (Pa$^2$)", fontsize=fontsize)
    contour(reshape(xobs, npts, npts), reshape(yobs, npts, npts), reshape(stressfield, npts, npts), stresscontours, linewidths = 0.5, colors = "k")
    for i in 1:els.endidx
        # if els.name[i] == "fault"
            plot([els.x1[i], els.x2[i]], [els.y1[i], els.y2[i]], "-k", linewidth = 2.0, zorder=40, color="lightgray")
        # end
    end
    xlim([minimum(xobs), maximum(xobs)])
    ylim([minimum(yobs), maximum(yobs)])
    xticks([minimum(xobs), 0, maximum(xobs)])
    yticks([minimum(yobs), 0, maximum(yobs)])
    gca().set_aspect("equal")
    xlabel(L"$x$ (m)", fontsize=fontsize)
    ylabel(L"$y$ (m)", fontsize=fontsize)
    ax.tick_params(labelsize = fontsize)
    suptitle(figuretitle, fontsize=fontsize)
    tight_layout()
    show()
    return nothing
end


function qdvis()
    foldername = "/home/meade/Desktop/data/qdvis/" * basename
    filenamejld2 =  foldername * "/" * basename * ".jld2"

    # filenamejld2 = basefilename * ".jld2"
    sol = load(filenamejld2, "integrator.sol")
    els = load(filenamejld2, "els")
    nu = load(filenamejld2, "nu")
    mu = load(filenamejld2, "mu")
    siay = 365.25 * 24 * 60 * 60

    # Create convenience structures
    idx = getidxdict(els)
    partialsconst = initpartials(els)

    # Create observation grid
    npts = 50
    xobs, yobs = obsgrid(-10e3, -5e3, 10e3, 5e3, npts)

    # Partials for BEM solve which I need for the surface
    @time _, _, partialsconst["trac"]["fault"]["fault"] = partialsconstdispstress(slip2dispstress, els, idx["fault"], idx["fault"], mu, nu)

    # Contours of base10 log of velocity magnitudes through time
    nfault = Integer(length(sol[1])[1] / 3)
    timeidx = 900
    # for i in timeidx:timeidx
    for i in 1:length(sol.t)
        @show i

        # Forward model for volumetric displacements and stresses
        dispfault, stressfault = constdispstress(slip2dispstress, xobs, yobs, els, idx["fault"], sol[i][1:3:end], sol[i][2:3:end], mu, nu)

        # Plots single time step
        figuretitle = @sprintf("time step = %06.0f, time = %08.2f", i, sol.t[i] / siay)
        close("all")
        plotlocal(els, idx, dispfault, stressfault, xobs, yobs, npts, figuretitle)

        # Save figure
        plt.savefig(@sprintf("%06.0f.png", i))
    end

end
qdvis()
# png2mp4(30, basename)
