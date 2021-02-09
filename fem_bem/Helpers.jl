using PyCall
using PyPlot
using FEniCS

function fenics_plot_scalar(f)
    py"""
    from fenics import *
    def FF(f):
        import matplotlib.pyplot as plt
        plt.figure()
        c = plot(f, colobar=True)
        plt.colorbar(c)
        plt.show()
    """
    py"FF"(f.pyobject)
end

function fenics_plot_vec2d(f)
    py"""
    from fenics import *
    def FF(f):
        import matplotlib.pyplot as plt
        fx = dot(f,Constant((1.0, 0.0)))
        fy = dot(f,Constant((0.0, 1.0)))
        plt.figure()
        plt.subplot(1,2,1)
        c = plot(fx, colobar=True)
        plt.colorbar(c)
        plt.subplot(1,2,2)
        c = plot(fy, colobar=True)
        plt.colorbar(c)
        plt.tight_layout()
        plt.show()
    """
    py"FF"(f.pyobject)
end

# Calculate strain
function fenics_strain(u)
    0.5 * (nabla_grad(u) + FEniCS.Transpose(nabla_grad(u)))
end


# Calculate
function fenics_stress(u, lambda, mu)
     lambda * nabla_div(u) * Identity(2) + 2 * mu * fenics_strain(u)
end

function fenics_solve(V, f)
    # Create mesh and define function space
 
    bc = DirichletBC(V, Constant((0, 0)), "on_boundary") # What BCs are bing set???

    # Solve
    u = TrialFunction(V)
    d = geometric_dimension(u) # space dimension
    v = TestFunction(V)
    T = Constant((0, 0))
    a = inner(fenics_stress(u, lambda, mu), fenics_strain(v)) * dx
    L = FEniCS.dot(f, v) * dx + FEniCS.dot(T, v) * ds
    u = FeFunction(V)
    lvsolve(a, L, u, bc)
    u
end

function fenics_interior_eval(f, x, y)
    py"""
    def G(f, x, y):
        return f(x, y)
    """
    py"G"(f.pyobject, x, y)
end

function fenics_tensor_components(V, f)
    fxx_fxy = project(FEniCS.dot(f, Constant((1,0))), V)
    fxy_fyy = project(FEniCS.dot(f, Constant((0,1))), V)
    fxx_fxy, fxy_fyy
end

function fenics_eval_u(fenics_u, idx, x, y)
    n = length(idx)
    u = zeros(2*n)
    for idx_i in 1:n
        i = idx[idx_i]
        u[(2*idx_i - 1):(2*idx_i)] = fenics_interior_eval(fenics_u, x[i], y[i])
    end
    u
end

function fenics_eval_s(sxx_sxy, sxy_syy, idx, x, y)
    n = length(idx)
    s = zeros(3*n)
    for idx_i in 1:n
        i = idx[idx_i]
        sxx, sxy = fenics_interior_eval(sxx_sxy, x[i], y[i])
        sxy2, syy = fenics_interior_eval(sxy_syy, x[i], y[i])
        s[3*idx_i - 2] = sxx
        s[3*idx_i - 1] = syy
        s[3*idx_i - 0] = sxy
    end
    s
end

function stress_to_trac(S, idx, xn, yn)
    n = length(idx)
    t = zeros(2*n)
    for idx_i in 1:n
        i = idx[idx_i]
        sxx, syy, sxy = S[3*idx_i-2:3*idx_i]
        tx = sxx * xn[i] + sxy * yn[i]
        ty = sxy * xn[i] + syy * yn[i]
        t[2*idx_i - 1] = tx
        t[2*idx_i] = ty
    end
    t
end

function twopanel(xgrid, ygrid, npts, U, S, idx, els; figsize=(12,3), ylim=[-20000, 2000])
    # Start of nice visualization
    figure(figsize=figsize)
    subplot(1, 2, 1)
    field = sqrt.(U[:, 1].^2 + U[:, 2].^2)
    ncontours = 10
    lowfield = 0
    highfield = 100
    ncontours = LinRange(lowfield, highfield, 11)    

    xlim = [-20000 20000]
    scale = 1.0
    fieldmax = maximum(@.abs(field))
    contourf(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
             reshape(field, npts, npts), ncontours,
             vmin=lowfield, vmax=highfield,
             cmap = get_cmap("magma"))
    clim(lowfield, highfield)
    colorbar(fraction=0.020, pad=0.05, extend = "both", label = L"$||u||$ (m)")
    contour(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
            reshape(field, npts, npts), ncontours,
            vmin=lowfield, vmax=highfield,
            linewidths=0.25, colors="w")
    plotelements(els)
    gca().set_aspect("equal")
    gca().set_xlim([xlim[1], xlim[2]])
    gca().set_ylim([ylim[1], ylim[2]])
    gca().set_xticks([-20000, -10000, 0, 10000, 20000])
    gca().set_yticks([-20000, -10000, 0])
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")
    xv = [els.x1[idx["T"]] ; els.x2[idx["T"]][end] ; 20000 ; -20000]
    yv = [els.y1[idx["T"]] ; els.y2[idx["T"]][end]; 2000 ; 2000]
    fill(xv, yv, "lightgray", zorder=30)

    # Find principle stress orientations
    stressdiff = zeros(length(U[:, 1]))
    for i in 1:length(U[:,1])
        mat = [S[i, 1] S[i, 3]; S[i, 3] S[i, 2]]
        ev = eigvals(mat)
        stressdiff[i] = log10(abs(ev[2]-ev[1]))
    end
    
    subplot(1, 2, 2)
    field = stressdiff
    xlim = [-20000 20000]
    lowfield = 5
    highfield = 10
    ncontours = LinRange(lowfield,highfield, 11)    
    contourf(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
             reshape(field, npts, npts), ncontours,
             vmin=lowfield, vmax=highfield,
             cmap=get_cmap("viridis"))
    clim(lowfield, highfield)
    colorbar(fraction=0.020, pad=0.05, extend="both",
             label=L"$\Delta \sigma$ (Pa)")
    contour(reshape(xgrid, npts, npts), reshape(ygrid, npts, npts),
            reshape(field, npts, npts), ncontours,
            vmin=lowfield, vmax=highfield,
            linewidths=0.25, colors="w")
    plotelements(els)
    gca().set_aspect("equal")
    gca().set_xlim([xlim[1], xlim[2]])
    gca().set_ylim([ylim[1], ylim[2]])
    gca().set_xticks([-20000, -10000, 0, 10000, 20000])
    gca().set_yticks([-20000, -10000, 0])
    xlabel(L"$x$ (m)")
    ylabel(L"$y$ (m)")
    xv = [els.x1[idx["T"]] ; els.x2[idx["T"]][end] ; 20000 ; -20000]
    yv = [els.y1[idx["T"]] ; els.y2[idx["T"]][end]; 2000 ; 2000]
    fill(xv, yv, "lightgray", zorder=30)
    tight_layout()
    return nothing
end