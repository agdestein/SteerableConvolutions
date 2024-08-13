# LSP hack                                     #src
if false                                       #src
    include("../src/SteerableConvolutions.jl") #src
    using .SteerableConvolutions               #src
end                                            #src

# # Scratch
#
# Scratchpad test script.

using SteerableConvolutions

G = CyclicGroup(8)
g = G(3)

r1 = Representation([0, 1, 2, 3, 4], randn(8, 8))
r2 = Representation([0, 1, 2], randn(5, 5))
r3 = irrep(G, 3)

directsum(map(i -> irrep(G, i), irreps(G))...)(elements(G)[1])
display.(map(i -> irrep(G, i)(elements(G)[1]), irreps(G)));

directsum(r1)
r = r1 ⊕ r2 ⊕ r3
r(elements(G)[1])

irrep(G, 1)

irrepmat(g, 0)

r1.basis
r1(g)

e = elements(G)
r = regular_representation(G)

r.irreps
r.basis

irrep(G, 4)(e[1])

r(e[5])
round.(Int, r(e[8]))

elements(G)
