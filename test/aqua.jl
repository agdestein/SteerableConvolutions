@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(SteerableConvolutions; ambiguities = false)
end
