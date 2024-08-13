# Add SteerableConvolutions environment,
# since test environment is not allowed to depend on the package directly
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using TestItemRunner

@run_package_tests
