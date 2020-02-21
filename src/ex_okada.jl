using Revise
using PyCall
using PyPlot
using Bem2d

function ex_okada()
    mu = 30e9
    nu = 0.25

    # Flat fault
    nfault = 1
    x1, y1, x2, y2 = Bem2d.discretizedline(-1, 0, 1, 0, nfault)
    els = Bem2d.Elements(Int(1e5))
    for i in 1:length(x1)
        els.x1[els.endidx + i] = x1[i]
        els.y1[els.endidx + i] = y1[i]
        els.x2[els.endidx + i] = x2[i]
        els.y2[els.endidx + i] = y2[i]
        els.name[els.endidx + i] = "fault"
    end
    Bem2d.standardize_elements!(els)
    idx = Bem2d.getidxdict(els)

    obswidth = 10
    npts = 100
    x, y = Bem2d.obsgrid(-obswidth, -obswidth, obswidth, obswidth, npts)

    # Constant slip fault
    strikeslip = 1.0 * ones(els.endidx)
    tensileslip = 0.0 * ones(els.endidx)
    UfU, σfU = Bem2d.constdispstress(slip2dispstress, x, y, els, idx["fault"], strikeslip, tensileslip, mu, nu)

    # Okada solution
    ow = pyimport("okada_wrapper") # from okada_wrapper import dc3dwrapper
    uokada = zeros(length(x), 2)
    σokada = zeros(length(x), 3)

    for i in 1:length(x)
        deep = 1000.0
        _, u, ∇u = ow.dc3dwrapper(0.6, [0.0, x[i], y[i]-deep], 0.0 + deep,
                                  0.0, [-0.7, 0.7], [-0.7, 0.7],
                                  [0.0, 0.0, 1.0])
        uokada[i, 1] = u[2]
        uokada[i, 2] = u[3]
    end
    PyPlot.close("all")
    Bem2d.plotfields(els, reshape(x, npts, npts), reshape(y, npts, npts),
               uokada, σokada, "Okada strike-slip")
    PyPlot.show()
end
ex_okada()
