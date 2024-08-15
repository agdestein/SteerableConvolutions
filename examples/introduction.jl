# LSP hack                                     #src
if false                                       #src
    include("../src/SteerableConvolutions.jl") #src
    using .SteerableConvolutions               #src
end                                            #src

# # Introduction
#
# Introduction to groups, representations, and steerable convolutions
# using the SteerableConvolutions.jl package.

using SteerableConvolutions

# ## Groups
#
# A group ``G = (E, *)`` is a set of elements ``E`` equipped with a group action ``*``.
# Consider the cyclic group ``G = C_8`` of order ``N = 8``.
# This group contains 8 elements ``E = \{0, 1, 2, 3, 4, 5, 6, 7\}``, and the
# group action ``*`` represents rotations in the plane of angle ``2 \pi n / N``
# for ``n \in E``.
G = CyclicGroup(8)

# We can access or create elements of the group.
# An element ``g`` is represented by an integer ``n`` between 0 and 7,
# representing a rotation by angle ``2 \pi n / N``.
g = G(3)

# A group element ``g`` has some properties:
typeof(g), g.group, g.n

# Since the group is finite, all the group elements can
# be accessed with the [`elements`](@ref) function:
elements(G)

# Group elements can be multiplied.
# This group action represents consecutive rotations, which creates a new group
# element.
# Note how the angles start over after 8 rotations.
G(1) * G(2), G(3) * G(7)

# The unit element ``e = 0`` represents a rotation by zero degrees:
one(G), one(g), G(0), G(2) * G(-2)

# The unit elements does not modify the angle:
one(G) * G(3)

# Since the group action is analogous to multiplication, we can compute inverses and powers:
inv(g), g^2

# ## Group representations
#
# A [`Representation`](@ref) of a group ``G`` is a matrix-valued function over the group:
#
# ```math
# \rho: G \to \mathbb{R}^{d \times d}
# ```
#
# where ``d`` is the size of the representation.
#
# ### Regular representation
#
# One such representation is the regular representation, which encodes
# the group action using a permutation matrix of size ``N``:
ρ = regular_representation(G)
round.(Int, ρ(g))

# Internally, all [`Representation`](@ref)s are stored in a different basis, leading
# to some round-off errors:
ρ(g)

# One property of group representations is that group actions are represented
# as matrix multiplications:
g, h = G(3), G(-1)
g * h

# The corresponding matrix is
matrix = ρ(g * h)
round.(Int, matrix)

# Compare with
matrix = ρ(g) * ρ(h)
round.(Int, matrix)

# As a result, the matrix of the unit element is always the identity matrix:
round.(Int, ρ(one(G)))

# ### Irreducible representations
#
# The regular representation is composed of [`IrreducibleRepresentation`](@ref)s (irreps).
# For the circular group, there are ``N / 2`` irreps of different frequencies.
# With frequency ``0``, we get the trivial (constant) representation ``\psi_0``,
# which is a ``1 \times 1`` identity matrix for all group elements:
ψ₀ = irrep(G, 0)
ψ₀(G(0)), ψ₀(G(1)), ψ₀(G(2))

# The next irrep ``\psi_1`` has frequency ``1`` and takes the value of a
# ``2 \times 2`` rotation matrix ``\psi_1(g) = R(2 \pi g / N)``, with
#
# ```math
# R(\theta) =
# \begin{bmatrix}
#     \cos(\theta) & -\sin(\theta) \\
#     \sin(\theta) & \cos(\theta)
# \end{bmatrix}.
# ```
ψ₁ = irrep(G, 1)
ψ₁(g)

# All the valid irrep frequencies are accessible through the `irreps` function:
irreps(G)

# ### Direct sums
#
# The direct sum of two represenations is a new representation that gives
# the block diagonal concatenation of the original matrices.
# It can be called using the [`directsum`](@ref) function or
# the [`⊕`](@ref) operator (type `\oplus<tab>`):
ρ_sum = ψ₀ ⊕ ψ₀ ⊕ ψ₀ ⊕ ψ₁
ρ_sum(G(1))

#-
ρ_sum(G(2))

#-
ρ_sum(G(3))

# Any representation ``\rho`` can be represented as a direct sum
# of irreps through a change of basis:
#
# ```math
# \rho(g) = B \left( \bigoplus_{i \in I} \psi_i(g) \right) B^{-1} \quad \forall g \in G,
# ```
#
# where ``B`` is the basis matrix and ``(\psi_i)_i`` are the irreps of frequencies ``i``.
# The list of indices ``I`` can contain repetitions.
# Let's check for the regular representation `ρ`:
B = basis(ρ)

# The irreps of ρ are stored as a list of frequencies ``I``:
I = irreps(ρ)

# We can assemble the direct sum:
ρ_block = ⊕(map(i -> irrep(G, i), I)...)
ρ_block(g)

# We can verify that the direct sum is equal to the regular representation
# when applying the change of basis:
B * ρ_block(g) * inv(B) ≈ ρ(g)
