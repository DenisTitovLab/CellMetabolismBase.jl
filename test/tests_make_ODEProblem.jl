@testitem "make_ODEProblem" setup=[TestMetabolicPathway] begin
    using LabelledArrays, BenchmarkTools, OrdinaryDiffEq

    test_pathway = TestMetabolicPathway.test_pathway

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

    function test_odes!(dmetabs, metabs, params, t)
        dmetabs.A_media = 0.0
        dmetabs.A = rate_enz1(metabs, params) - rate_enz2(metabs, params)
        dmetabs.B = 2 * rate_enz2(metabs, params) - rate_enz3(metabs, params)
        dmetabs.C = rate_enz3(metabs, params) - 2 * rate_enz4(metabs, params)
        dmetabs.D = rate_enz4(metabs, params)
        return nothing
    end

    metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0, D = 1.0)
    dmetabs = similar(metabs)
    params = LVector(
        Enz1_Vmax = 1.0,
        Enz1_Keq = 10.0,
        Enz1_K_A_media = 1.0,
        Enz1_K_A = 1.0,
        Enz2_Vmax = 1.0,
        Enz2_Keq = 1.0,
        Enz2_K_A = 2.0,
        Enz2_K_B = 1.0,
        Enz3_Vmax = 1.0,
        Enz3_Keq = 10.0,
        Enz3_K_B = 1.0,
        Enz3_K_C = 1.0,
        Enz4_Vmax = 1.0,
        Enz4_Keq = 1.0,
        Enz4_K_C = 1.0,
        Enz4_K_D = 1.0,
    )

    dmetabs_test_calc = @LArray rand(5) (:A_media, :A, :B, :C, :D)
    rates = (1.0, 0.5, 0.3, 0.2)  # Rates for Enz1, Enz2, Enz3, Enz4
    CellMetabolismBase.calculate_dmetabs_from_enz_rates!(
        test_pathway,
        dmetabs_test_calc,
        rates,
    )
    @test dmetabs_test_calc.A_media ≈ 0.0  # Constant metabolite should not change
    @test dmetabs_test_calc.A ≈ 1.0 - 0.5  # Produced by Enz1, consumed by Enz2
    @test dmetabs_test_calc.B ≈ 2 * 0.5 - 0.3  # Produced 2x by Enz2, consumed by Enz3
    @test dmetabs_test_calc.C ≈ 0.3 - 2 * 0.2  # Produced by Enz3, consumed 2x by Enz4
    @test dmetabs_test_calc.D ≈ 0.2  # Produced by Enz4

    dmetabs = @LArray rand(5) (:A_media, :A, :B, :C, :D)
    dmetabs_expected = @LArray rand(5) (:A_media, :A, :B, :C, :D)
    dmetabs_wrong_names = @LArray rand(5) (:A_media, :A, :B, :C, :X)
    @test_throws Exception CellMetabolismBase.metabolicpathway_odes!(test_pathway, dmetabs_wrong_names, metabs, params, 0.0)
    CellMetabolismBase.metabolicpathway_odes!(test_pathway, dmetabs, metabs, params, 0.0)
    test_odes!(dmetabs_expected, metabs, params, 0.0)
    @test dmetabs == dmetabs_expected

    benchmark_result = @benchmark CellMetabolismBase.metabolicpathway_odes!(
        $test_pathway,
        $dmetabs,
        $metabs,
        $params,
        0.0,
    )
    @test mean(benchmark_result.times) <= 200 #ns
    @test benchmark_result.allocs == 0

    benchmark_result = @benchmark test_odes!($dmetabs, $metabs, $params, 0.0)
    @test mean(benchmark_result.times) <= 200 #ns
    @test benchmark_result.allocs == 0

    prob_manual = ODEProblem(test_odes!, metabs, (0.0, 1e6), params)
    prob = make_ODEProblem(test_pathway, metabs, (0.0, 1e6), params)
    @test prob isa ODEProblem
    @test prob.tspan == (0.0, 1e6)
    @test prob.u0 == metabs
    @test prob.p == params

    enzymes = CellMetabolismBase._generate_Enzymes(test_pathway)
    manual_benchmark_result = @benchmark CellMetabolismBase._generate_Enzymes($test_pathway)
    @test mean(manual_benchmark_result.times) <= 10 #ns
    @test manual_benchmark_result.allocs == 0
    @test enzymes isa Tuple{Vararg{Enzyme}}
    @test length(enzymes) == 4
    @test enzymes[1] == Enzyme(:Enz1,(:A_media,),(:A,))
    @test enzymes[2] == Enzyme(:Enz2,(:A,),(:B, :B))
    @test enzymes[3] == Enzyme(:Enz3,(:B,),(:C,))
    @test enzymes[4] == Enzyme(:Enz4,(:C, :C),(:D,))
    mismatched_metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0)
    @test_throws Exception make_ODEProblem(metab_pathway, mismatched_metabs, tspan, params)

    sol_manual = OrdinaryDiffEq.solve(
        prob_manual,
        Rodas5P(),
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )
    sol = OrdinaryDiffEq.solve(
        prob,
        Rodas5P(),
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )
    @test sol_manual == sol

    sol_manual = OrdinaryDiffEq.solve(
        prob_manual,
        OrdinaryDiffEq.Rodas4P(),
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )
    sol = OrdinaryDiffEq.solve(
        prob,
        OrdinaryDiffEq.Rodas4P(),
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )
    @test sol_manual == sol

    manual_benchmark_result = @benchmark OrdinaryDiffEq.solve(
        prob_manual,
        Rodas5P(),
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )
    @test mean(manual_benchmark_result.times) <= 20_000_000 #ns
    @test manual_benchmark_result.allocs < 30_000

    package_benchmark_result = @benchmark OrdinaryDiffEq.solve(
        prob,
        Rodas5P(),
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )
    @test mean(package_benchmark_result.times) <= 20_000_000 #ns
    @test package_benchmark_result.allocs < 30_000

    @test manual_benchmark_result.allocs >= package_benchmark_result.allocs

    println(
        "Manual ODE solve=$(round(mean(manual_benchmark_result.times)/1_000_000; sigdigits=2))ms",
    )
    println("Manual ODE allocations=$(manual_benchmark_result.allocs)")
    println(
        "Package ODE solve=$(round(mean(package_benchmark_result.times)/1_000_000; sigdigits=2))ms",
    )
    println("Package ODE allocations=$(package_benchmark_result.allocs)")
end
