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
field = FiberField(fieldtype, x)

# We now apply a group transform (rotation by 90 degrees):
g = G(1)
newfield = g * field

# Plot the all fields
function plot(field)
    (; x) = field
    s = @. sqrt(field.x[:, :, 2]^2 + field.x[:, :, 3]^2)
    color = s[:] ./ maximum(s)
    limits = (0, 1, 0, 1)
    fig = Figure(; size = (800, 400))
    heatmap(fig[1, 1], grid, grid, x[:, :, 1]; axis = (; limits, title = "Scalar field"))
    arrows(
        fig[1, 2],
        grid,
        grid,
        x[:, :, 2],
        x[:, :, 3];
        lengthscale = 0.1,
        color,
        axis = (; limits, title = "Vector field"),
    )
    fig
end

# Here are the original fields:
plot(field)

# Here are the rotated fields. For the vector field, note also how both the domain and the arrows have been rotated.
plot(newfield)

# ## Equivariance
#
# A function ``f`` mapping one fiber field ``u`` with representation ``\rho_\text{in}``
# to another fiber field of representation ``\rho_\text{out}`` is
# said to be ``G``-equivariant if
#
# ```math
#     f(\rho_\text{in}(g) u) = \rho_\text{out}(g) f(u)
# ```
#
# for all ``g \in G``.
#
# For example, consider ``\rho_\text{in}`` to be the representation from above,
# and ``\rho_\text{out}`` to be the trivial representation:
fieldtype_out = FieldType(gspace, [irrep(gspace.group, 0)])

# The following function is not equivariant:
function notequivariant(u)
    (; x) = u
    y = x[:, :, 2:2] + x[:, :, 3:3]
    FiberField(fieldtype_out, y)
end

# since
a = g * notequivariant(field)
b = notequivariant(g * field)
a.x ≈ b.x

# is `false`. The following function is equivariant:
function norm2(u)
    (; x) = u
    y = @. x[:, :, 2:2]^2 + x[:, :, 3:3]^2
    FiberField(fieldtype_out, y)
end

# since
a = g * norm2(field)
b = norm2(g * field)
a.x ≈ b.x
