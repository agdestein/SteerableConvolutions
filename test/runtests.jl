push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Aqua
using SteerableConvolutions
using Test

@testset "SteerableConvolutions" begin
    include("aqua.jl")
end
