using Revise
using AbstractPlotting
using Makie

npts = 100
limits = FRect(-50, 0, 200, 1)
x = collect(LinRange(0, 10, npts))
yupdate = Node(rand(npts))

scene = plot(x, yupdate, limits=limits)
plot!(scene, x, rand(npts), color = :green, limits=limits)
display(scene)

for i=1:100
    yupdate[] = rand(npts)
    sleep(1/50)
end
