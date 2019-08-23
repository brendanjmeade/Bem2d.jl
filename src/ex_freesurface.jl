using Revise
using PyCall
using PyPlot
using Bem2d

function ex_freesurface()
    # Material properties and observation grid
    μ = 30e9
    ν = 0.25
    npts = 50
    obswidth = 5
    xobs, yobs = obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Free surface
    elements = Elements()
    x1, y1, x2, y2 = discretizedline(-5, 0, 5, 0, 20)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "freesurface"
    end
    standardize_elements!(elements)

    # 45 degree dipping fault
    x1, y1, x2, y2 = discretizedline(-1, -1, 0, 0, 1)
    for i in 1:length(x1)
        elements.x1[i + elements.endidx] = x1[i]
        elements.y1[i + elements.endidx] = y1[i]
        elements.x2[i + elements.endidx] = x2[i]
        elements.y2[i + elements.endidx] = y2[i]
        elements.name[i + elements.endidx] = "fault"
    end
    standardize_elements!(elements)

    # Constant slip fault
    srcidx = findall(x->x == "fault", elements.name)
    obsidx = findall(x->x == "freesurface", elements.name)
    u1, σ1, t1 = ∂constslip(elements, srcidx, obsidx, μ, ν)
    srcidx = findall(x->x == "freesurface", elements.name)
    obsidx = findall(x->x == "freesurface", elements.name)
    u2, σ2, t2 = ∂constslip(elements, srcidx, obsidx, μ, ν)

    # Constant case: Predict surface displacements from unit strike slip forcing
    faultslip = [1 ; 1]
    ufullspace = u1 * faultslip
    ufreesurface = inv(t2) * (t1 * faultslip)
    plot(ufullspace[1:2:end], "r.")
    plot(ufreesurface[1:2:end], "b+")
    
    # gca().set_aspect("equal")
    # gca().set_xlim([xlim[1], xlim[2]])
    # gca().set_ylim([ylim[1], ylim[2]])
    # gca().set_xticks([xlim[1], xlim[2]])
    # gca().set_yticks([ylim[1], ylim[2]])
    xlabel("x (m)")
    ylabel("u (m)")
    show()

    # plot(dispfreesurface[1:2:end])
    # plot!(dispfreesurface[2:2:end])
end
ex_freesurface()