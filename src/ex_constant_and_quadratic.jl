using Revise
using Bem2d

# Create a flat fault
elements = Elements();
L = 1e3;
x1, y1, x2, y2 = discretizedline(-L, 0, L, 0, 3)
for i in 1:length(x1)
    push!(elements.x1, x1[i]);
    push!(elements.y1, y1[i]);
    push!(elements.x2, x2[i]);
    push!(elements.y2, y2[i]);
    push!(elements.name, "fault01");
end

standardize_elements!(elements)
