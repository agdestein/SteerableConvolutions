@testitem "CyclicGroup" begin
    G = CyclicGroup(8)
    @test G(8) == G(0)
    @test G(-1) == G(7)
    @test one(G) * G(3) == G(3)
    @test G(1) * G(2) == G(3)
    @test inv(G(3)) * G(3) == one(G)
end

@testitem "Representation" begin
    G = CyclicGroup(8)
    I = [0, 0, 1, 4]
    B = randn(5, 5)
    ρ = Representation(I, B)
    @test frequencies(ρ) == I
    @test basis(ρ) == B
    g, h = G(3), G(2)
    @test ρ(one(G)) ≈ one(ρ(g))
    @test ρ(inv(g)) ≈ inv(ρ(g))
    @test ρ(g * h) ≈ ρ(g) * ρ(h)
end

@testitem "Irrep frequencies" begin
    @test frequencies(CyclicGroup(8)) == 0:4
    @test frequencies(CyclicGroup(9)) == 0:4
    G = CyclicGroup(8)
    @test all(==([1;;]), irrep(G, 0).(elements(G)))
end

@testitem "Direct sum" begin
    G = CyclicGroup(8)
    ρ1 = Representation([0, 0, 1, 4], randn(5, 5))
    ρ2 = Representation([4], randn(1, 1))
    ρ3 = irrep(G, 1)
    ρ = ρ1 ⊕ ρ2 ⊕ ρ3
    @test frequencies(ρ) == vcat(frequencies(ρ1), frequencies(ρ2), frequencies(ρ3))
end
