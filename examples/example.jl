if false
    include("../src/SteerableConvolutions.jl")
    using .SteerableConvolutions
end

# # Example
#
# Example run.

using SteerableConvolutions

group = CyclicGroup(8)
one(group)

r1 = Representation([0, 1, 2, 3, 4], randn(8, 8))
r2 = Representation([0, 1, 2], randn(5, 5))
r3 = irrep(group, 3)

directsum(map(i -> irrep(group, i), irreps(group))...)(elements(group)[1])
display.(map(i -> irrep(group, i)(elements(group)[1]), irreps(group)));

directsum(r1)
r = r1 ⊕ r2 ⊕ r3
r(elements(group)[1])

irrep(group, 1)

g = group(3)

group(2) * group(3)

irrepmat(g, 0)

r1.basis
r1(g)

e = elements(group)
r = regular_representation(group)


r.irreps
r.basis

irrep(group, 4)(e[1])

r(e[5])
round.(Int, r(e[8]))

elements(group)
