# LSP hack                                     #src
if false                                       #src
    include("../src/SteerableConvolutions.jl") #src
    using .SteerableConvolutions               #src
end                                            #src

# # Equivariance
#
# Introduction to [`GSpace`](@ref)s, [`FieldType`](@ref)s, and [`FiberField`](@ref)s
# using the SteerableConvolutions.jl package.

using SteerableConvolutions
using CairoMakie

# ## GSpaces
#
# Create 2-dimensional [`GSpace`](@ref) with cyclic group of order 4.
G = CyclicGroup(4)
gspace = GSpace(G, 2)

# ## Field types
#
# Create a [`FieldType`](@ref) from a list of representations.
representations = [irrep(G, 0), irrep(G, 1)]
fieldtype = FieldType(gspace, representations)

# ## Fiber fields
#
# Since our fieldtype contains one scalar field (from `irrep(G, 0)`) and
# one 2D-vector field (from `irrep(G, 1)`), we need to create a feature
# field with three channels.
# We wrap it in a [`FiberField`](@ref) from a field type and data.
n = 50
grid = LinRange(0, 1, n + 1)
mask = @. (1 - grid^2) * (1 - (grid')^2)
x = zeros(n + 1, n + 1, 3)
@. x[:, :, 1] = sinpi(10 * grid) * mask
@. x[:, :, 2] = sinpi(2 * grid) * cospi(3 * grid') * mask
@. x[:, :, 3] = -cospi(2 * grid) * sinpi(3 * grid') * mask
f = FiberField(fieldtype, x)

# We now apply a group transform (rotation by 90 degrees):
g = G(1)
fnew = g * f

# Plot the all fields
function plot(f)
    s = @. sqrt(f.x[:, :, 2]^2 + f.x[:, :, 3]^2)
    color = s[:] ./ maximum(s)
    fig = Figure(; size = (800, 400))
    heatmap(
        fig[1, 1],
        grid,
        grid,
        f.x[:, :, 1];
        axis = (; limits = (0, 1, 0, 1), title = "Scalar field"),
    )
    arrows(
        fig[1, 2],
        grid,
        grid,
        f.x[:, :, 2],
        f.x[:, :, 3];
        lengthscale = 0.1,
        color,
        axis = (; limits = (0, 1, 0, 1), title = "Vector field"),
    )
    fig
end

# Here are the original fields:
plot(f)

# Here are the rotated fields. For the vector field, note also how both the domain and the arrows have been rotated.
plot(fnew)
