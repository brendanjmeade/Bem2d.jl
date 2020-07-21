using Revise
using PyPlot
using Infiltrator
using Bem2d


"""
    dislocationinabox()

Comparing half-space and dislocaiton in a box solutions
"""
function dislocationinabox()
    close("all")
    mu = 30e9
    nu = 0.25
    alpha = 1e-7 # scalar preconditioner

    # Element geometries and data structures for half space approximation
    els = Elements(Int(1e5))
    nfault = 20
    x1, y1, x2, y2 = discretizedline(-10e3, -10e3, 0, 0, nfault) # 45 degree dipping fault
    addelsez!(els, x1, y1, x2, y2, "fault")
    nsurf = 50
    x1, y1, x2, y2 = discretizedline(-50e3, 0, 50e3, 0, nsurf) # Free surface
    addelsez!(els, x1, y1, x2, y2, "surf")
    idx = getidxdict(els)

    # Solve the BEM problem
    _, H_surf_fault = PUTC(slip2dispstress, els, idx["surf"], idx["fault"], mu, nu)
    _, H_surf_surf = PUTC(slip2dispstress, els, idx["surf"], idx["surf"], mu, nu)
    faultslip = sqrt(2) / 2 * ones(2 * nfault)
    @infiltrate
    Ueff = inv(H_surf_surf) * (H_surf_fault * faultslip)

    # Forward evaluation
    Ufault, _ = constdispstress(
        slip2dispstress,
        els.xcenter[idx["surf"]],
        els.ycenter[idx["surf"]],
        els,
        idx["fault"],
        ones(nfault),
        ones(nfault),
        mu,
        nu,
    )
    Usurf, _ = constdispstress(
        slip2dispstress,
        els.xcenter[idx["surf"]],
        els.ycenter[idx["surf"]],
        els,
        idx["surf"],
        Ueff[1:2:end],
        Ueff[2:2:end],
        mu,
        nu,
    )
    Utotal = @. Ufault - Usurf

    # Alternative BEM solution
    T_fault_all, _ =
        PUTC(slip2dispstress, els, idx["fault"], collect(1:1:els.endidx), mu, nu)
    _, H_surf_all = PUTC(slip2dispstress, els, idx["surf"], collect(1:1:els.endidx), mu, nu)
    bcs = zeros(2 * els.endidx)
    bcs[1:2*nfault] .= -1
    Ueffalt = inv([T_fault_all; H_surf_all]) * bcs

    # Alternative forward evaluation
    Ualt, _ = constdispstress(
        slip2dispstress,
        els.xcenter[idx["surf"]],
        els.ycenter[idx["surf"]],
        els,
        collect(1:1:els.endidx),
        Ueffalt[1:2:end],
        Ueffalt[2:2:end],
        mu,
        nu,
    )

    # Plot surface displacements
    fontsize = 20
    markersize = 12
    linewidth = 2.0
    figure(figsize = (12, 12))
    ax = subplot(2, 1, 1) # x-displacements
    plot(els.xcenter[idx["surf"]], Ueff[1:2:end], "b.", label = "Ueff")
    plot(els.xcenter[idx["surf"]], Ueffalt[2*nfault+1:2:end], "g.", label = "Ueffalt")
    plot(els.xcenter[idx["surf"]], Ualt[1:2:end], "c.", label = "Ualt")
    plot(els.xcenter[idx["surf"]], Utotal[:, 1], "r.", label = "total")
    legend()
    xlabel(L"$x$ (m)")
    ylabel(L"$u_x$ (m)")

    ax = subplot(2, 1, 2) # y-displacements
    plot(els.xcenter[idx["surf"]], Ueff[2:2:end], "b.", label = "Ueff")
    plot(els.xcenter[idx["surf"]], Ueffalt[2*nfault+2:2:end], "g.", label = "Ueffalt")
    plot(els.xcenter[idx["surf"]], Ualt[2:2:end], "c.", label = "Ualt")
    plot(els.xcenter[idx["surf"]], Utotal[:, 2], "r.", label = "total")
    legend()
    xlabel(L"$x$ (m)")
    ylabel(L"$u_y$ (m)")
end
dislocationinabox()
