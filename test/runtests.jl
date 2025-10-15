# using CellMetabolismBase
# using Test
# using SafeTestsets
# @safetestset "Test make_ODEProblem" begin
#     include("tests_make_ODEProblem.jl")
# end

include("support/enzyme_rate_fixtures.jl")

using TestItemRunner
@run_package_tests
