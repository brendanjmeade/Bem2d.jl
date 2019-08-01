using Revise
using Bem2d

# Create a flat fault
elements = Elements()
L = 1e3
x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, 3)
for i in 1:length(x1)
    elements.x1[i + elements.lastidx] = x1[i]
    elements.y1[i + elements.lastidx] = y1[i]
    elements.x2[i + elements.lastidx] = x2[i]
    elements.y2[i + elements.lastidx] = y2[i]
    elements.uxconst[i + elements.lastidx] = 1.0
    elements.uyconst[i + elements.lastidx] = 1.0
    elements.uxquad[i + elements.lastidx, :] = [1.0 1.0 1.0]
    elements.uyquad[i + elements.lastidx, :] = [1.0 1.0 1.0]    
    elements.name[i + elements.lastidx] = "fault01"
end
standardize_elements!(elements)
