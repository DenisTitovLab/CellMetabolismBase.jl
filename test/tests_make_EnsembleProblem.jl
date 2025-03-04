@testitem "make_EnsembleProblem" begin
    using LabelledArrays, BenchmarkTools, OrdinaryDiffEq

    # Use the same metabolic pathway as in tests_make_ODEProblem.jl
    test_pathway = MetabolicPathway(
        (:A_media,),
        (
            (:Enz1, (:A_media,), (:A,)),
            (:Enz2, (:A,), (:B, :B,)),
            (:Enz3, (:B,), (:C,)),
            (:Enz4, (:C, :C), (:D,)),
        ),
    )

    # Define rate functions for enzymes
    function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz1,(:A_media,),(:A,)}, metabs, params)
        return params.Enz1_Vmax * (metabs.A_media - metabs.A / params.Enz1_Keq) /
               (1 + metabs.A_media / params.Enz1_K_A_media + metabs.A / params.Enz1_K_A)
    end

    function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz2,(:A,),(:B, :B)}, metabs, params)
        return params.Enz2_Vmax * (metabs.A - metabs.B^2 / params.Enz2_Keq) /
               (1 + metabs.A / params.Enz2_K_A + metabs.B^2 / params.Enz2_K_B)
    end

    function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz3,(:B,),(:C,)}, metabs, params)
        return params.Enz3_Vmax * (metabs.B - metabs.C / params.Enz3_Keq) /
               (1 + metabs.B / params.Enz3_K_B + metabs.C / params.Enz3_K_C)
    end

    function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz4,(:C, :C),(:D,)}, metabs, params)
        return params.Enz4_Vmax * (metabs.C^2 - metabs.D / params.Enz4_Keq) /
               (1 + metabs.C^2 / params.Enz4_K_C + metabs.D / params.Enz4_K_D)
    end

    # Create test initial conditions and parameters
    init_cond1 = LVector(A_media=2.0, A=1.0, B=1.0, C=1.0, D=1.0)
    init_cond2 = LVector(A_media=3.0, A=0.5, B=0.5, C=0.5, D=0.5)

    params1 = LVector(
        Enz1_Vmax=1.0, Enz1_Keq=10.0, Enz1_K_A_media=1.0, Enz1_K_A=1.0,
        Enz2_Vmax=1.0, Enz2_Keq=1.0, Enz2_K_A=2.0, Enz2_K_B=1.0,
        Enz3_Vmax=1.0, Enz3_Keq=10.0, Enz3_K_B=1.0, Enz3_K_C=1.0,
        Enz4_Vmax=1.0, Enz4_Keq=1.0, Enz4_K_C=1.0, Enz4_K_D=1.0
    )

    params2 = LVector(
        Enz1_Vmax=2.0, Enz1_Keq=8.0, Enz1_K_A_media=0.5, Enz1_K_A=0.5,
        Enz2_Vmax=1.5, Enz2_Keq=2.0, Enz2_K_A=1.0, Enz2_K_B=2.0,
        Enz3_Vmax=0.8, Enz3_Keq=12.0, Enz3_K_B=1.5, Enz3_K_C=0.8,
        Enz4_Vmax=1.2, Enz4_Keq=0.5, Enz4_K_C=0.8, Enz4_K_D=1.2
    )

    # Test case: Main function with multiple initial conditions and parameters
    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        [init_cond1, init_cond2],
        [params1, params2]
    )

    @test ensemble_prob isa EnsembleProblem

    # Test solving the ensemble problem
    sol = solve(ensemble_prob, RadauIIA9(), EnsembleSerial(); abstol=1e-15, reltol=1e-8, trajectories=2)
    @test length(sol) == 2

    # Check that the solutions are different (used different params/init conditions)
    @test sol[1].u[end] != sol[2].u[end]

    # Test case: Single parameter set with multiple initial conditions
    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        [init_cond1, init_cond2],
        params1
    )

    @test ensemble_prob isa EnsembleProblem

    sol = solve(ensemble_prob, RadauIIA9(), EnsembleSerial(); abstol=1e-15, reltol=1e-8, trajectories=2)
    @test length(sol) == 2

    # Solutions should differ due to different initial conditions
    @test sol[1].u[1] != sol[2].u[1]

    # Test case: Single initial condition with multiple parameter sets
    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        init_cond1,
        [params1, params2]
    )

    @test ensemble_prob isa EnsembleProblem

    sol = solve(ensemble_prob, RadauIIA9(), EnsembleSerial(); abstol=1e-15, reltol=1e-8, trajectories=2)
    @test length(sol) == 2

    # Solutions should differ due to different parameters
    @test sol[1].u[end] != sol[2].u[end]

    # Test case: Custom time span
    tspan = (0.0, 100.0)
    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        [init_cond1, init_cond2],
        [params1, params2],
        tspan=tspan
    )

    @test ensemble_prob.prob.tspan == tspan

    # Test case: Error handling - mismatched vector lengths
    @test_throws ErrorException make_EnsembleProblem(
        test_pathway,
        [init_cond1, init_cond2],
        [params1]
    )
end
