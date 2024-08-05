# # Example
#
# Example run.

using SteerableConvolutions

G = CyclicGroup{8}
r1 = Representation([0, 1, 2, 3, 4], randn(8, 8))
r2 = Representation([0, 1, 2], randn(5, 5))
r3 = irrep(G, 3)

directsum(map(i -> irrep(G, i), irreps(G))...)(elements(G)[8])

directsum(r1)
r = r1 ⊕ r2 ⊕ r3
r(elements(G)[1])

G = CyclicGroup{8}
irrep(G, 0)

g = G(3)

irrepmat(g, 0)

r1.basis
r1(g)

G = CyclicGroup{8}
e = elements(G)
r = regular_representation(G)

r.irreps
r.basis

irrep(G, 4)(e[7])

round.(Int, r(e[5]))

elements(G)
