# Set up examples environment

using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(; path = joinpath(@__DIR__, "..")))
Pkg.instantiate()
