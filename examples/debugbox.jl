using Revise
using PyPlot
using Infiltrator
using LinearAlgebra
using Bem2d


"""
    plotparticular(g, rho, lambda, mu)

Plot fundamental particular fields
"""
function plotparticular(g, rho, lambda, mu)
    npts = 100
    x, y = obsgrid(-20000, -20000, 20000, 20000, npts) 
    Up, Sp = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    figure(figsize=(20, 10))
    subplot(2, 3, 1)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Up[:, 1], npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("ux")
    subplot(2, 3, 2)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Up[:, 2], npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("uy")
    subplot(2, 3, 4)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Sp[:, 1], npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("Sxx")
    subplot(2, 3, 5)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Sp[:, 2], npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("Syy")
    subplot(2, 3, 6)
    contourf(reshape(x, npts, npts), reshape(y, npts, npts), reshape(Sp[:, 3], npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("Sxy")
end


"""
    gravitysquareparticular()

Experiments with gravity body force.
"""
function gravitysquareparticular()
    # TODO: It's strange that the top of the model has to be at zero.  Can we generalize this?
    # TODO: Rule of thumb for choosing precondtioner value (alpha)?
    # TODO: Make a version of gravityparticularfunctions() for BC generation

    close("all")

    alpha = 7e-8 # scalar preconditioner for traction terms
    fontsize = 20
    mu = 3e10
    lambda = 3e10
    nu = 0.25
    rho = 2700
    g = 9.81
    nels = 20
    npts = 100
    L = 1e4
    offset = 10
    # x, y = obsgrid(-L+offset, -2*L+offset, L-offset, 0-offset, npts) 

    plotparticular(g, rho, lambda, mu)
    return

    # Define BEM geometry
    els = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-L, -2*L, L, -2*L, nels) # Bottom
    addelsez!(els, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(L, -2*L, L, 0, nels) # Right hand side
    addelsez!(els, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(L, 0, -L, 0, nels) # Top
    addelsez!(els, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-L, 0, -L, -2*L, nels) # Left hand side
    addelsez!(els, x1, y1, x2, y2, "L")

    # Common indexing
    idx = getidxdict(els) # Should this return "all" - YES TODO
    bcidxU = idx["B"] # Boundaries with *displacement* BCs
    bcidxT = [idx["R"] ; idx["T"]; idx["L"]] # Boundaries with *traction* BCs
    bcidxall = collect(1:1:els.endidx) # All boundaries

    # Gravity square problem with quadratic elements
    T_pU_qall, _ = PUTQ(slip2dispstress, els, bcidxU, bcidxall, mu, nu)
    _, H_pT_qall = PUTQ(slip2dispstress, els, bcidxT, bcidxall, mu, nu)
    TH = [T_pU_qall ; alpha .* H_pT_qall] # Assemble combined linear operator

    # Particular solution and effective boundary conditions
    xnodes = transpose(els.xnodes[idx["B"], :])[:]
    ynodes = transpose(els.ynodes[idx["B"], :])[:]
    UB, _ = gravityparticularfunctions(xnodes, ynodes, g, rho, lambda, mu)
    bcs = zeros(6 * els.endidx)
    bcs[1:2:6*nels] = UB[:, 1] # Bottom boundary (x-component)
    bcs[2:2:6*nels] = UB[:, 2] # Bottom boundary (y-component)
    bcs *= -1 # This is neccesary for the right answer and is consistent with derivation

    # BEM solve to get particular solution
    Ueffparticular = inv(TH) * bcs

    # Evaluate and plot interior solution
    Uinteriorcomplementary, Sinteriorcomplementary = quaddispstress(slip2dispstress, x, y, els, bcidxall, quadstack(Ueffparticular[1:2:end]), quadstack(Ueffparticular[2:2:end]), mu, nu)
    Uinteriorparticular, Sinteriorparticular = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    U = @. Uinteriorcomplementary + Uinteriorparticular
    S = @. Sinteriorcomplementary + Sinteriorparticular
    Umag = sqrt.(U[:, 1].^2 + U[:, 2].^2)

    figure()
    contourf(reshape(x, npts, npts), reshape(y, npts, npts),
             reshape(Umag, npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("BEM - below ground")

    #####
    ##### Now with top of box not at zero
    #####
    # New observaiton coordinates
    x, y = obsgrid(-L+offset, 0, L-offset, 2*L-offset, npts)

    # Define BEM geometry
    # elsbox = Elements(Int(1e5))
    # x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, nels) # Bottom
    # addelsez!(elsbox, x1, y1, x2, y2, "B")
    # x1, y1, x2, y2 = discretizedline(L, 0, L, 2*L, nels) # Right hand side
    # addelsez!(elsbox, x1, y1, x2, y2, "R")
    # x1, y1, x2, y2 = discretizedline(L, 2*L, -L, 2*L, nels) # Top
    # addelsez!(elsbox, x1, y1, x2, y2, "T")
    # x1, y1, x2, y2 = discretizedline(-L, 2*L, -L, 0, nels) # Left hand side
    # addelsez!(elsbox, x1, y1, x2, y2, "L")
    # Define BEM geometry
    num = -2*L
    y = y .+ num
    elsbox = Elements(Int(1e5))
    x1, y1, x2, y2 = discretizedline(-L, 0+num, L, 0+num, nels) # Bottom
    addelsez!(elsbox, x1, y1, x2, y2, "B")
    x1, y1, x2, y2 = discretizedline(L, 0+num, L, 2*L+num, nels) # Right hand side
    addelsez!(elsbox, x1, y1, x2, y2, "R")
    x1, y1, x2, y2 = discretizedline(L, 2*L+num, -L, 2*L+num, nels) # Top
    addelsez!(elsbox, x1, y1, x2, y2, "T")
    x1, y1, x2, y2 = discretizedline(-L, 2*L+num, -L, 0+num, nels) # Left hand side
    addelsez!(elsbox, x1, y1, x2, y2, "L")

    # Common indexing
    idxbox = getidxdict(elsbox) # Should this return "all" - YES TODO
    bcidxU = idxbox["B"] # Boundaries with *displacement* BCs
    bcidxT = [idxbox["R"] ; idxbox["T"]; idxbox["L"]] # Boundaries with *traction* BCs
    bcidxall = collect(1:1:elsbox.endidx) # All boundaries

    T_B_BRTL, H_B_BRTL = PUTC(slip2dispstress, elsbox, idxbox["B"],
                              [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                              mu, nu)
    T_R_BRTL, H_R_BRTL = PUTC(slip2dispstress, elsbox, idxbox["R"],
                              [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                              mu, nu)
    T_L_BRTL, H_L_BRTL = PUTC(slip2dispstress, elsbox, idxbox["L"],
                              [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                              mu, nu)
    T_T_BRTL, H_T_BRTL = PUTC(slip2dispstress, elsbox, idxbox["T"],
                              [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                              mu, nu)
    THbox = [T_B_BRTL ; H_R_BRTL; H_T_BRTL ; H_L_BRTL]
    
    # Displacement BCs for bottom
    Ugabove, Sgabove = gravityparticularfunctions(elsbox.xcenter[idxbox["B"]],
                                                  elsbox.ycenter[idxbox["B"]],
                                                  g, rho, lambda, mu)
    bcsboxgravity = zeros(2 * elsbox.endidx)
    bcsboxgravity[1:2:2*nels] = Ugabove[:, 1]
    bcsboxgravity[2:2:2*nels] = Ugabove[:, 2]

    # Traction BCs for top and sides
    # Ug, Sg = gravityparticularfunctions(elsbox.xcenter[idxbox["T"]],
    #                                     elsbox.ycenter[idxbox["T"]],
    #                                     g, rho, lambda, mu)  
    # Ttx = zeros(nels)
    # Tty = zeros(nels)
    # for i in 1:length(Sg[:, 1])
    #     nvec = [elsbox.xnormal[idxbox["T"][i]] ; elsbox.ynormal[idxbox["T"][i]]]
    #     temp = [Sg[i, 1] Sg[i, 3] ; Sg[i, 3] Sg[i, 2]] * nvec
    #     Ttx[i] = temp[1]
    #     Tty[i] = temp[2]
    # end
    # bcsboxgravity[4*nels+1:2:6*nels] = Ttx
    # bcsboxgravity[4*nels+2:2:6*nels] = Tty
    bcsboxgravity *= -1
    Ueffboxparticular = THbox \ bcsboxgravity
    
    # Forward solution on grid
    Ucomp, Scomp = constdispstress(slip2dispstress, x, y, elsbox,
                                   [idxbox["B"] ; idxbox["R"] ; idxbox["T"] ; idxbox["L"]],
                                   Ueffboxparticular[1:2:end], Ueffboxparticular[2:2:end], mu, nu)
    Uint, Sint = gravityparticularfunctions(x, y, g, rho, lambda, mu)
    Ugravityonly = @. Ucomp + Uint
    Sgravityonly = @. Scomp + Sint
    bemaboveumag = sqrt.(Ugravityonly[:, 1].^2 + Ugravityonly[:, 2].^2)

    figure();
    contourf(reshape(x, npts, npts), reshape(y, npts, npts),
             reshape(bemaboveumag, npts, npts), 50)
    colorbar()
    gca().set_aspect("equal")
    title("BEM - above ground")
end
gravitysquareparticular()
