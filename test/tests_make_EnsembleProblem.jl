@testitem "make_EnsembleProblem" setup=[TestMetabolicPathway] begin
    using LabelledArrays, BenchmarkTools, OrdinaryDiffEqFIRK, SciMLBase, Distributions

    # Use the shared metabolic pathway
    test_pathway = TestMetabolicPathway.test_pathway

    # Create test initial conditions and parameters
    init_cond1 = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0, D = 1.0)
    init_cond2 = LVector(A_media = 3.0, A = 0.5, B = 0.5, C = 0.5, D = 0.5)

    params1 = LVector(
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

    params2 = LVector(
        Enz1_Vmax = 2.0,
        Enz1_Keq = 8.0,
        Enz1_K_A_media = 0.5,
        Enz1_K_A = 0.5,
        Enz2_Vmax = 1.5,
        Enz2_Keq = 2.0,
        Enz2_K_A = 1.0,
        Enz2_K_B = 2.0,
        Enz3_Vmax = 0.8,
        Enz3_Keq = 12.0,
        Enz3_K_B = 1.5,
        Enz3_K_C = 0.8,
        Enz4_Vmax = 1.2,
        Enz4_Keq = 0.5,
        Enz4_K_C = 0.8,
        Enz4_K_D = 1.2,
    )

    # Test case: Main function with multiple initial conditions and parameters
    ensemble_prob =
        make_EnsembleProblem(test_pathway, [init_cond1, init_cond2], [params1, params2])

    @test ensemble_prob isa EnsembleProblem

    # Test solving the ensemble problem
    sol = solve(
        ensemble_prob,
        RadauIIA9(),
        EnsembleSerial();
        abstol = 1e-15,
        reltol = 1e-8,
        trajectories = 2,
    )
    @test length(sol) == 2

    # Check that the solutions are different (used different params/init conditions)
    @test sol.u[1].u[end] != sol.u[2].u[end]

    # Test case: Single parameter set with multiple initial conditions
    ensemble_prob = make_EnsembleProblem(test_pathway, [init_cond1, init_cond2], params1)

    @test ensemble_prob isa EnsembleProblem

    sol = solve(
        ensemble_prob,
        RadauIIA9(),
        EnsembleSerial();
        abstol = 1e-15,
        reltol = 1e-8,
        trajectories = 2,
    )
    @test length(sol) == 2

    # Solutions should differ due to different initial conditions
    @test sol.u[1].u[1] != sol.u[2].u[1]

    # Test case: Single initial condition with multiple parameter sets
    ensemble_prob = make_EnsembleProblem(test_pathway, init_cond1, [params1, params2])

    @test ensemble_prob isa EnsembleProblem

    sol = solve(
        ensemble_prob,
        RadauIIA9(),
        EnsembleSerial();
        abstol = 1e-15,
        reltol = 1e-8,
        trajectories = 2,
    )
    @test length(sol) == 2

    # Solutions should differ due to different parameters
    @test sol.u[1].u[end] != sol.u[2].u[end]

    # Test case: Custom time span
    tspan = (0.0, 100.0)
    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        [init_cond1, init_cond2],
        [params1, params2],
        tspan = tspan,
    )

    @test ensemble_prob.prob.tspan == tspan

    # Test case: Error handling - mismatched vector lengths
    @test_throws ErrorException make_EnsembleProblem(
        test_pathway,
        [init_cond1, init_cond2],
        [params1],
    )

    # Create distributions for initial conditions
    init_cond_dist = LVector(
        A_media = Uniform(1.5, 2.5),
        A = Normal(1.0, 0.1),
        B = Normal(1.0, 0.1),
        C = Normal(1.0, 0.1),
        D = Normal(1.0, 0.1),
    )

    # Create distributions for parameters
    params_dist = LVector(
        Enz1_Vmax = Uniform(0.8, 1.2),
        Enz1_Keq = Normal(10.0, 1.0),
        Enz1_K_A_media = LogNormal(0.0, 0.2),
        Enz1_K_A = LogNormal(0.0, 0.2),
        Enz2_Vmax = Uniform(0.8, 1.2),
        Enz2_Keq = Normal(1.0, 0.2),
        Enz2_K_A = LogNormal(2.0, 0.2),
        Enz2_K_B = LogNormal(0.0, 0.2),
        Enz3_Vmax = Uniform(0.8, 1.2),
        Enz3_Keq = Normal(10.0, 1.0),
        Enz3_K_B = LogNormal(0.0, 0.2),
        Enz3_K_C = LogNormal(0.0, 0.2),
        Enz4_Vmax = Uniform(0.8, 1.2),
        Enz4_Keq = Normal(1.0, 0.2),
        Enz4_K_C = LogNormal(0.0, 0.2),
        Enz4_K_D = LogNormal(0.0, 0.2),
    )

    # Test with a small number of bootstraps for faster testing
    n_bootstraps = 100
    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        init_cond_dist,
        params_dist;
        n_bootstraps = n_bootstraps,
    )

    # Basic type check
    @test ensemble_prob isa EnsembleProblem

    # Test solving the ensemble problem
    sol = solve(
        ensemble_prob,
        RadauIIA9(),
        EnsembleSerial();
        abstol = 1e-12,
        reltol = 1e-8,
        trajectories = n_bootstraps,
    )
    @test length(sol) == n_bootstraps

    # Check that the solutions vary (due to different random samples)
    final_vals = [s.u[end][4] for s in sol.u]  # Get final values for metabolite C
    @test length(unique(final_vals)) > 1   # Should have different values

    # Test with custom time span
    tspan = (0.0, 100.0)
    ensemble_prob2 = make_EnsembleProblem(
        test_pathway,
        init_cond_dist,
        params_dist;
        n_bootstraps = n_bootstraps,
        tspan = tspan,
    )
    @test ensemble_prob2.prob.tspan == tspan

    # Test that random sampling works properly by checking that
    # sampled values are within expected ranges
    prob_func = ensemble_prob.prob_func
    base_prob = ensemble_prob.prob

    # Check a few samples
    for i = 1:5
        new_prob = prob_func(base_prob, i, 1)

        # Test initial conditions are within reasonable bounds
        # For normal distributions, most values should be within 3 std devs
        @test 1.5 <= new_prob.u0[1] <= 2.5  # A_media ~ Uniform(1.5, 2.5)
        @test 0.7 <= new_prob.u0[2] <= 1.3  # A ~ Normal(1.0, 0.1)

        # Test parameters are within reasonable bounds
        @test 0.8 <= new_prob.p[1] <= 1.2   # Enz1_Vmax ~ Uniform(0.8, 1.2)
        @test 7.0 <= new_prob.p[2] <= 13.0  # Enz1_Keq ~ Normal(10.0, 1.0)
    end

    # Test case: Initial conditions as distributions with fixed parameters
    init_cond_dist = LVector(
        A_media = Uniform(1.5, 2.5),
        A = Normal(1.0, 0.1),
        B = Normal(1.0, 0.1),
        C = Normal(1.0, 0.1),
        D = Normal(1.0, 0.1),
    )

    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        init_cond_dist,
        params1;
        n_bootstraps = n_bootstraps,
    )

    @test ensemble_prob isa EnsembleProblem

    # Test solving the ensemble problem
    sol = solve(
        ensemble_prob,
        RadauIIA9(),
        EnsembleSerial();
        abstol = 1e-12,
        reltol = 1e-8,
        trajectories = n_bootstraps,
    )
    @test length(sol) == n_bootstraps

    # Check that the solutions vary due to different initial conditions
    initial_vals = [s.u[1][2] for s in sol.u]  # Get initial values for metabolite A
    @test length(unique(initial_vals)) > 1  # Should have different values

    # Verify that all trajectories used the same parameters
    prob_func = ensemble_prob.prob_func
    base_prob = ensemble_prob.prob
    for i = 1:5
        new_prob = prob_func(base_prob, i, 1)
        @test new_prob.p == params1  # Parameters should remain fixed
    end

    # Test case: Fixed initial conditions with parameters as distributions
    ensemble_prob = make_EnsembleProblem(
        test_pathway,
        init_cond1,
        params_dist;
        n_bootstraps = n_bootstraps,
    )

    @test ensemble_prob isa EnsembleProblem

    # Test solving the ensemble problem
    sol = solve(
        ensemble_prob,
        RadauIIA9(),
        EnsembleSerial();
        abstol = 1e-12,
        reltol = 1e-8,
        trajectories = n_bootstraps,
    )
    @test length(sol) == n_bootstraps

    # Check that the solutions vary due to different parameters
    final_vals = [s.u[end][5] for s in sol.u]  # Get final values for metabolite D
    @test length(unique(final_vals)) > 1  # Should have different values

    # Verify that all trajectories used the same initial conditions
    prob_func = ensemble_prob.prob_func
    base_prob = ensemble_prob.prob
    for i = 1:5
        new_prob = prob_func(base_prob, i, 1)
        @test new_prob.u0 == init_cond1  # Initial conditions should remain fixed
    end
end
