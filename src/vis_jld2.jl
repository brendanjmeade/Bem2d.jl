using Revise
using JLD2
using PyCall
using PyPlot

PyPlot.close("all")

function plottimeseries()
    filename = "2019-10-18T12:02:03.37.jld2"
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
    filename = "2019-10-18T12:02:03.37.jld2"
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

    # Calculate cumulative slip
    sx = zeros(size(vx))
    # for itime in 1:length(t)
    #     for ifault in 1:nfault
    #         # sx[itime, ifault] = sum(vx[1:itime, ifault])
    #         sx[itime, ifault] = sum(vx[1:itime, ifault])
    #     end
    # end

    # Swap order of loops here?
    # dt = diff(sol.t)
    for ifault in 1:nfault
        for itime in 2:length(t)
            # sx[itime, ifault] = sx[itime - 1, ifault] + vx[itime, ifault] * dt[itime - 1]
            sx[itime, ifault] = sx[itime - 1, ifault] + vx[itime, ifault] * (t[itime] - t[itime - 1])
        end
    end


    # PyPlot.figure(figsize = (15, 5))
    # plotme = log10.(vx')
    # PyPlot.contourf(plotme, 200, cmap = get_cmap("plasma"))
    # PyPlot.colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}v$ (m/s)")
    # PyPlot.contour(plotme, 20, linewidths = 0.5, linestyles = "solid", colors = "k")
    # PyPlot.xlabel("time step")
    # PyPlot.ylabel("element index")
    # PyPlot.show()

    PyPlot.figure(figsize = (15, 5))
    # for i in 1:30
    #     PyPlot.plot(sx[], "-k")
    # end
    plotme = sx'
    PyPlot.contourf(plotme, 200, cmap = get_cmap("plasma"))
    PyPlot.colorbar(fraction = 0.020, pad = 0.05, extend = "both", label = L"$\log_{10}v$ (m/s)")
    PyPlot.contour(plotme, 20, linewidths = 0.5, linestyles = "solid", colors = "k")
    PyPlot.xlabel("time step")
    PyPlot.ylabel("element index")
    PyPlot.show()

end
plotvelocities()
