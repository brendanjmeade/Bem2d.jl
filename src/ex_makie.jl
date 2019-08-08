using Makie

y = range(-0.997669, stop = 0.997669, length = 23)
contour(range(-0.99, stop = 0.99, length = 23), y, rand(23, 23), levels = 10)
