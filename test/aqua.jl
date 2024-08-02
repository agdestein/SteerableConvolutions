@testset "Aqua" begin
    @info "Testing code with Aqua"
    Aqua.test_all(SteerableConvolutions; ambiguities = false)
end
