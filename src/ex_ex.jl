using Revise
using Makie

# # TODO: Add lines instead of updating
# x = collect(LinRange(0, 4, 400))
# y = rand(400)
# scene = lines(x, y)
# lineplot = scene[end]
# display(scene)
# for i in 1:50
#     y = rand(400)
#     lines!(scene, x, y, color=:green)
#         # lineplot[2] = y
#     AbstractPlotting.update!(scene)
#     sleep(0.01)
# end


x = collect(LinRange(0, 1, 100))
scene = lines(x, rand(length(x)), color = :blue)
display(scene)

for i in 1:10
    lines!(scene, x, rand(length(x)), color = :red)
    display(scene)
    # AbstractPlotting.update!(scene)
end
