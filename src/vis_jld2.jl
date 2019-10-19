using Revise
using JLD2
using PyCall
using PyPlot
using Bem2d

PyPlot.close("all")

function plottimeseries()
    filename = "2019-10-18T14:40:57.735.jld2"
    JLD2.@load filename sol

    siay = 365.25 * 24 * 60 * 60
    nfault = Int(size(sol.u[1])[1] / 3)
    t = [x / siay for x in sol.t]
    θ = zeros(length(t), nfault)
    vx = zeros(length(t), nfault)
    vy = zeros(length(t), nfault)
    for i in 1:length(t)
        vx[i, :] = sol.u[i][1:3:end]
        vy[i, :] = sol.u[i][2:3:end]
        θ[i, :] = sol.u[i][3:3:end]
    end

    PyPlot.figure(figsize = (15, 8))

    PyPlot.subplot(3, 2, 1)
    PyPlot.plot(t, vx, "-", linewidth = 0.5)
    PyPlot.yscale("log")
    PyPlot.ylabel(L"v_x")
    PyPlot.subplot(3, 2, 3)
    PyPlot.plot(t, vy, "-", linewidth = 0.5)
    PyPlot.yscale("log")
    PyPlot.ylabel(L"v_y")
    PyPlot.subplot(3, 2, 5)
    PyPlot.plot(t, θ, "-", linewidth = 0.5)
    PyPlot.yscale("log")
    PyPlot.xlabel("t (years)")
    PyPlot.ylabel(L"\theta")

    PyPlot.subplot(3, 2, 2)
    PyPlot.plot(1:1:length(t), vx, "-", linewidth = 0.5)
    PyPlot.yscale("log")
    PyPlot.ylabel(L"v_x")
    PyPlot.subplot(3, 2, 4)
    PyPlot.plot(1:1:length(t), vy, "-", linewidth = 0.5)
    PyPlot.yscale("log")
    PyPlot.ylabel(L"v_y")
    PyPlot.subplot(3, 2, 6)
    PyPlot.plot(1:1:length(t), θ, "-", linewidth = 0.5)
    PyPlot.yscale("log")
    PyPlot.xlabel("time step #")
    PyPlot.ylabel(L"\theta")

    PyPlot.figure(figsize = (15, 5))
    plotme = log10.(vx')
    PyPlot.contourf(plotme, 200, cmap = get_cmap("plasma"))
    PyPlot.colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}v$ (m/s)")
    PyPlot.contour(plotme, 20, linewidths = 0.5, linestyles = "solid", colors = "k")
    PyPlot.xlabel("time step")
    PyPlot.ylabel("element index")
    PyPlot.show()
end
# plottimeseries()

function plotvelocities()
    filename = "2019-10-18T14:40:57.735.jld2"
    # filename = "2019-10-18T12:02:03.37.jld2"

    JLD2.@load filename sol

    blockvelx = 1e-9 # should read this in from .jld2 file
    siay = 365.25 * 24 * 60 * 60
    nfault = Int(size(sol.u[1])[1] / 3)
    t = [x / siay for x in sol.t]
    θ = zeros(length(t), nfault)
    vx = zeros(length(t), nfault)
    vy = zeros(length(t), nfault)
    for i in 1:length(t)
        vx[i, :] = sol.u[i][1:3:end]
        vy[i, :] = sol.u[i][2:3:end]
        θ[i, :] = sol.u[i][3:3:end]
    end

    # Calculate cumulative slip
    sx = zeros(size(vx))
    sxlinear = zeros(size(vx))

    for itime in 2:length(t)
        for ifault in 1:nfault
            sx[itime, ifault] = sx[itime - 1, ifault] + vx[itime, ifault] * (t[itime] - t[itime - 1])
            sxlinear[itime, ifault] = sxlinear[itime - 1, ifault] + blockvelx * (t[itime] - t[itime - 1])
        end
    end

    fontsize = 24
    PyPlot.figure(figsize = (35, 7))
    plotme = log10.(transpose(vx))
    PyPlot.contourf(plotme, 8, cmap = rycroftcmap())
    cbar = PyPlot.colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}v$ (m/s)")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.ax.set_ylabel(L"$\log_{10}v_x$ (m/s)", fontsize = fontsize)
    PyPlot.contour(plotme, 8, linewidths = 0.5, linestyles = "solid", colors = "k")
    PyPlot.xlabel("time step", fontsize=fontsize)
    PyPlot.ylabel("element index", fontsize=fontsize)
    PyPlot.gca().tick_params(labelsize = fontsize)

    PyPlot.figure(figsize = (35, 7))
    plotme = transpose(sx)
    PyPlot.contourf(plotme, 20, cmap = get_cmap("plasma"))
    cbar = PyPlot.colorbar(fraction = 0.020, pad = 0.05, extend = "both")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.ax.set_ylabel(L"$s_x$ (m)", fontsize = fontsize)
    PyPlot.contour(plotme, 20, linewidths = 0.5, linestyles = "solid", colors = "k")
    PyPlot.xlabel("time step", fontsize=fontsize)
    PyPlot.ylabel("element index", fontsize=fontsize)
    PyPlot.gca().tick_params(labelsize = fontsize)

    PyPlot.figure(figsize = (35, 7))
    plotme = transpose(sxlinear - sx)
    PyPlot.contourf(plotme, 5, cmap = get_cmap("PRGn"), vmin=-1e-7, vmax=1e-7)
    cbar = PyPlot.colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$block motion - slip$ (m/s)")
    cbar.ax.tick_params(labelsize = fontsize)
    cbar.ax.set_ylabel(L"$\Delta s_x$ (m)", fontsize = fontsize)
    PyPlot.contour(plotme, 5, linewidths = 0.5, linestyles = "solid", colors = "k")
    PyPlot.xlabel("time step", fontsize=fontsize)
    PyPlot.ylabel("element index", fontsize=fontsize)
    PyPlot.gca().tick_params(labelsize = fontsize)

    PyPlot.show()

end
plotvelocities()
