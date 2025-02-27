using CellMetabolismBase
using Test
using Aqua
using JET

@testset "CellMetabolismBase.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CellMetabolismBase)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(CellMetabolismBase; target_defined_modules = true)
    end
    # Write your tests here.
end
