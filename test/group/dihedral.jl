@testitem "Elements" begin
    G = DihedralGroup(8)
    for f in (false, true)
        @test G(f, 8) == G(f, 0)
        @test G(f, 8) == G(f, 0)
        @test G(f, -1) == G(f, 7)
        @test one(G) * G(f, 3) == G(f, 3)
        @test inv(G(f, 3)) * G(f, 3) == one(G)
    end
    @test G(false, 1) * G(false, 2) == G(false, 3)
    @test G(false, 1) * G(true, 2) == G(true, 3)
    @test G(true, 1) * G(false, 2) == G(true, 7)
    @test G(true, 1) * G(true, 2) == G(false, 7)
end

@testitem "Representation" begin
    G = DihedralGroup(8)
    I = frequencies(G)
    n = sum(i -> size(Irrep(G, i)(one(G)), 1), I)
    B = randn(n, n)
    ρ = Representation(I, B)
    @test frequencies(ρ) == I
    @test basis(ρ) == B
    e = rand(elements(G), 5)
    for g in e
        @test ρ(one(G)) ≈ one(ρ(g))
        @test ρ(inv(g)) ≈ inv(ρ(g))
        for h in e
            @test ρ(g * h) ≈ ρ(g) * ρ(h)
        end
    end
end

@testitem "Irrep frequencies" begin
    @test frequencies(DihedralGroup(8)) ==
          [(false, 0), (true, 0), (true, 1), (true, 2), (true, 3), (true, 4), (false, 4)]
    @test frequencies(DihedralGroup(9)) ==
          [(false, 0), (true, 0), (true, 1), (true, 2), (true, 3), (true, 4)]

    # Test trivial representation
    G = DihedralGroup(8)
    @test all(==([1;;]), Irrep(G, (false, 0)).(elements(G)))
end

@testitem "Direct sum" begin
    G = DihedralGroup(8)
    ρ1 = Representation([(false, 0), (true, 0), (true, 1), (false, 4)], randn(5, 5))
    ρ2 = Representation([(true, 4)], randn(1, 1))
    ρ3 = Irrep(G, (true, 1))
    ρ = ρ1 ⊕ ρ2 ⊕ ρ3
    @test frequencies(ρ) == vcat(frequencies(ρ1), frequencies(ρ2), frequencies(ρ3))
end
