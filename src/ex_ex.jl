using Revise
using Makie

#     Makie.contour(range(-0.99, stop = 0.99, length = 23), y, rand(23, 23), levels = 10)

# function makeheatmaps(bufs, rows, columns)
#     heatmaps = []
#     n = rows*columns
#     @assert length(bufs) == n
#     for buf in bufs
#         # [end] necessary because heatmap returns the containing scene, not the heatmap
#         # hm = heatmap(buf, padding=(0,0), colorrange=(-60,-10), alpha=1.0)
#         hm = heatmap(buf, padding=(0,0), colorrange=(0,1), alpha=1.0)
#         axis = hm[Axis]
#         axis[:names][:axisnames] = ("", "")
#         push!(heatmaps, hm)
#     end
#
#     scene = AbstractPlotting.hbox(
#         [AbstractPlotting.vbox(heatmaps[(1:columns).+(row-1)*columns]...)
#          for row in 1:rows]...)
#     scene, [hm[end] for hm in heatmaps]
# end
#
# function runupdate(nframes)
#     datarows = 500
#     datacols = 500
#     plotrows = 4
#     plotcols = 4
#     n = plotrows*plotcols
#     bufs = [fill(0.0f0, datarows, datacols) for _ in 1:n]
#     scene, hms = makeheatmaps(bufs, plotrows, plotcols)
#     display(scene)
#
#     for _ in 1:nframes
#         frametime = @elapsed begin
#             for (hm, buf) in zip(hms, bufs)
#                 buf .= rand.(Float32)
#                 hm[:color] = buf
#             end
#         end
#         @show frametime
#         yield()
#     end
# end
# runupdate(100)



function makeheatmaps(bufs)
    rows = 1
    columns = 1
    heatmaps = []
    for buf in bufs
        # [end] necessary because heatmap returns the containing scene, not the heatmap
        # hm = heatmap(buf, padding=(0,0), colorrange=(-60,-10), alpha=1.0)
        hm = heatmap(buf, padding=(0,0), colorrange=(0,1), alpha=1.0)
        axis = hm[Axis]
        push!(heatmaps, hm)
    end

    scene = AbstractPlotting.hbox(
        [AbstractPlotting.vbox(heatmaps[(1:columns).+(row-1)*columns]...)
         for row in 1:rows]...)
    scene, [hm[end] for hm in heatmaps]
end

function runupdate(nframes)
    datarows = 500
    datacols = 500
    bufs = [fill(0.0f0, datarows, datacols)]
    scene, hms = makeheatmaps(bufs)

    display(scene)
    for _ in 1:nframes
        for (hm, buf) in zip(hms, bufs)
            buf .= rand.(Float32)
            hm[:color] = buf
        end
        yield()
    end
end
runupdate(100)
