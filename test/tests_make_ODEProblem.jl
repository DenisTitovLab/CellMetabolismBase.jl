@testitem "make_ODEProblem" begin
    using LabelledArrays, BenchmarkTools, OrdinaryDiffEq

    test_pathway = MetabolicPathway{
        (:A_media,),
        (
            (:Enz1, (:A_media,), (:A,)),
            (:Enz2, (:A,), (:B,:B,)),
            (:Enz3, (:B,), (:C,)),
            (:Enz4, (:C,:C), (:D,)),
        ),
    }()

    function CellMetabolismBase.rate(enzyme::Enzyme{:Enz1,(:A_media,),(:A,)}, metabs, params)
        return params.Enz1_Vmax * (metabs.A_media - metabs.A / params.Enz1_Keq) /
               (1 + metabs.A_media / params.Enz1_K_A_media + metabs.A / params.Enz1_K_A)
    end
    function CellMetabolismBase.rate(enzyme::Enzyme{:Enz2,(:A,),(:B,:B)}, metabs, params)
        return params.Enz2_Vmax * (metabs.A - metabs.B^2 / params.Enz2_Keq) /
               (1 + metabs.A / params.Enz2_K_A + metabs.B^2 / params.Enz2_K_B)
    end
    function CellMetabolismBase.rate(enzyme::Enzyme{:Enz3,(:B,),(:C,)}, metabs, params)
        return params.Enz3_Vmax * (metabs.B - metabs.C / params.Enz3_Keq) /
               (1 + metabs.B / params.Enz3_K_B + metabs.C / params.Enz3_K_C)
    end
    function CellMetabolismBase.rate(enzyme::Enzyme{:Enz4,(:C,:C),(:D,)}, metabs, params)
        return params.Enz4_Vmax * (metabs.C^2 - metabs.D / params.Enz4_Keq) /
               (1 + metabs.C^2 / params.Enz4_K_C + metabs.D / params.Enz4_K_D)
    end


    function rate_enz1(metabs, params)
        return params.Enz1_Vmax * (metabs.A_media - metabs.A / params.Enz1_Keq) /
               (1 + metabs.A_media / params.Enz1_K_A_media + metabs.A / params.Enz1_K_A)
    end
    function rate_enz2(metabs, params)
        return params.Enz2_Vmax * (metabs.A - metabs.B^2 / params.Enz2_Keq) /
               (1 + metabs.A / params.Enz2_K_A + metabs.B^2 / params.Enz2_K_B)
    end
    function rate_enz3(metabs, params)
        return params.Enz3_Vmax * (metabs.B - metabs.C / params.Enz3_Keq) /
               (1 + metabs.B / params.Enz3_K_B + metabs.C / params.Enz3_K_C)
    end
    function rate_enz4(metabs, params)
        return params.Enz4_Vmax * (metabs.C^2 - metabs.D / params.Enz4_Keq) /
               (1 + metabs.C^2 / params.Enz4_K_C + metabs.D / params.Enz4_K_D)
    end

    function test_odes(dmetabs, metabs, params, t)
        dmetabs.A_media = 0.0
        dmetabs.A = rate_enz1(metabs, params) - rate_enz2(metabs, params)
        dmetabs.B = 2*rate_enz2(metabs, params) - rate_enz3(metabs, params)
        dmetabs.C = rate_enz3(metabs, params) - 2*rate_enz4(metabs, params)
        dmetabs.D = rate_enz4(metabs, params)
        return nothing
    end

    metabs = LVector(A_media=2.0, A=1.0, B=1.0, C=1.0, D=1.0)
    dmetabs = similar(metabs)
    params = LVector(
        Enz1_Vmax=1.0,
        Enz1_Keq=10.0,
        Enz1_K_A_media=1.0,
        Enz1_K_A=1.0,
        Enz2_Vmax=1.0,
        Enz2_Keq=1.0,
        Enz2_K_A=2.0,
        Enz2_K_B=1.0,
        Enz3_Vmax=1.0,
        Enz3_Keq=10.0,
        Enz3_K_B=1.0,
        Enz3_K_C=1.0,
        Enz4_Vmax=1.0,
        Enz4_Keq=1.0,
        Enz4_K_C=1.0,
        Enz4_K_D=1.0,
    )

    benchmark_result = @benchmark CellMetabolismBase.metabolicpathway_odes!($test_pathway, $dmetabs, $metabs, $params, 0.0)
    @test mean(benchmark_result.times) <= 200 #ns
    @test benchmark_result.allocs == 0

    benchmark_result = @benchmark test_odes($dmetabs, $metabs, $params, 0.0)
    @test mean(benchmark_result.times) <= 200 #ns
    @test benchmark_result.allocs == 0

    prob_manual = ODEProblem(test_odes, metabs, (0.0, 1e6), params)
    prob = make_ODEProblem(test_pathway, metabs, (0.0, 1e6), params)
    @test prob isa ODEProblem

    sol_manual = solve(prob_manual, RadauIIA9(), abstol=1e-15, reltol=1e-8, save_everystep=false)
    sol = solve(prob, RadauIIA9(), abstol=1e-15, reltol=1e-8, save_everystep=false)
    @test sol_manual == sol

    sol_manual = solve(prob_manual, Rodas5P(), abstol=1e-15, reltol=1e-8, save_everystep=false)
    sol = solve(prob, Rodas5P(), abstol=1e-15, reltol=1e-8, save_everystep=false)
    @test sol_manual == sol

    manual_benchmark_result = @benchmark OrdinaryDiffEq.solve(prob_manual, RadauIIA9(), abstol=1e-15, reltol=1e-8, save_everystep=false)
    @test mean(manual_benchmark_result.times) <= 20_000_000 #ns
    @test manual_benchmark_result.allocs < 20_000

    package_benchmark_result = @benchmark OrdinaryDiffEq.solve(prob, RadauIIA9(), abstol=1e-15, reltol=1e-8, save_everystep=false)
    @test mean(package_benchmark_result.times) <= 20_000_000 #ns
    @test package_benchmark_result.allocs < 20_000

    @test manual_benchmark_result.allocs >= package_benchmark_result.allocs

    println("Manual ODE solve=$(round(mean(manual_benchmark_result.times)/1_000_000; sigdigits=2))ms")
    println("Manual ODE allocations=$(manual_benchmark_result.allocs)")
    println("Package ODE solve=$(round(mean(package_benchmark_result.times)/1_000_000; sigdigits=2))ms")
    println("Package ODE allocations=$(package_benchmark_result.allocs)")

end
