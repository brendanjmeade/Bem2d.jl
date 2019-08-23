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
    els = Elements()
    nfreesurf = 200
    x1, y1, x2, y2 = discretizedline(-5, 0, 5, 0, nfreesurf)
    els.x1[els.endidx + 1:els.endidx + nfreesurf] = x1
    els.y1[els.endidx + 1:els.endidx + nfreesurf] = y1
    els.x2[els.endidx + 1:els.endidx + nfreesurf] = x2
    els.y2[els.endidx + 1:els.endidx + nfreesurf] = y2
    els.name[els.endidx + 1:els.endidx + nfreesurf] .= "freesurf"
    standardize_elements!(els)

    # 45 degree dipping fault
    nfault = 1
    x1, y1, x2, y2 = discretizedline(-1, -1, 0, 0, nfault)
    els.x1[els.endidx + 1:els.endidx + nfault] = x1
    els.y1[els.endidx + 1:els.endidx + nfault] = y1
    els.x2[els.endidx + 1:els.endidx + nfault] = x2
    els.y2[els.endidx + 1:els.endidx + nfault] = y2
    els.name[els.endidx + 1:els.endidx + nfault] .= "fault"
    standardize_elements!(els)

    # Constant slip fault
    srcidx = findall(x->x == "fault", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u1, ∂σ1, ∂t1 = ∂constslip(els, srcidx, obsidx, μ, ν)
    srcidx = findall(x->x == "freesurf", els.name)
    obsidx = findall(x->x == "freesurf", els.name)
    ∂u2, ∂σ2, ∂t2 = ∂constslip(els, srcidx, obsidx, μ, ν)

    # Constant case: Predict surface displacements from unit strike slip forcing
    faultslip = sqrt(2) / 2 * [1 ; 1]
    ufullspace = ∂u1 * faultslip
    ufreesurface = inv(∂t2) * (∂t1 * faultslip)
    xplotconst = els.xcenter[findall(x->x == "freesurf", els.name)]
  
    close("all")
    figure(figsize = (8, 10))
    subplot(3, 1, 1)
    plotelements(els)
    gca().set_xlim([-5, 5])
    gca().set_xticks([-5, 0, 5])
    gca().set_aspect("equal")
    ylabel("y (m)")
    title("geometry")
  
    subplot(3, 1, 2)
    plot(xplotconst, ufullspace[1:2:end], "r-", label = "full space")
    plot(xplotconst, ufreesurface[1:2:end], "b-", label = "half space")    
    gca().set_xlim([-5, 5])
    gca().set_ylim([-0.6, 0.6])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([-0.5, 0.0, 0.5])
    legend()
    ylabel("u (m)")
    title("x displacements")

    subplot(3, 1, 3)
    plot(xplotconst, ufullspace[2:2:end], "r-", label = "full space")
    plot(xplotconst, ufreesurface[2:2:end], "b-", label = "half space")    
    gca().set_xlim([-5, 5])
    gca().set_ylim([-0.6, 0.6])
    gca().set_xticks([-5, 0, 5])
    gca().set_yticks([-0.5, 0.0, 0.5])
    legend()
    xlabel("x (m)")
    ylabel("u (m)")
    title("y displacements")
    show()
end
ex_freesurface()