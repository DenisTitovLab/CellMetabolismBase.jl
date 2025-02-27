using CellMetabolismBase
using Test
using TestItemRunner

@testitem "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(CellMetabolismBase)
end
@testitem "Code linting (JET.jl)" begin
    using JET
    JET.test_package(CellMetabolismBase; target_defined_modules=true)
end

@run_package_tests
