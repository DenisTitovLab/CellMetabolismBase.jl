@testitem "Accessor functions for MetabolicPathways" setup = [TestMetabolicPathway] begin
    using BenchmarkTools, LabelledArrays

    # Test basic Enzyme constructor and its errors
    @test Enzyme(:Enz1, (:A,), (:B,)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,)) === Enzyme{:Enz1,(:A,),(:B,),(),()}()
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1,), (:Inh1,)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1, :Act2), (:Inh1, :Inh2)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (), (:Inh1,)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1,), ()) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1,), (:Inh1,)) ===
          Enzyme{:Enz1,(:A,),(:B,),(:Act1,),(:Inh1,)}()
    @test Enzyme(:Enz1, (:A,), (:B,), (), ()) === Enzyme(:Enz1, (:A,), (:B,))
    @test_throws ErrorException Enzyme(:Enz1, (:A,), (:B,), [:Act1], (:Inh1,))
    @test_throws ErrorException Enzyme(:Enz1, (:A,), (:B,), (:Act1,), [:Inh1])
    @test_throws ErrorException Enzyme(:Enz1, (:A,), (:B,), "Act1", (:Inh1,))
    @test_throws ErrorException Enzyme(:Enz1, (:A,), (:B,), (:Act1,), "Inh1")
    @test_throws ErrorException Enzyme(:Enz1, [:A], (:B, :B))
    @test_throws ErrorException Enzyme(:Enz1, (:A,), [:B, :B])
    @test_throws ErrorException Enzyme((:Enz1,), (:A,), (:B,))

    # Test MetabolicPathway constructor and its errors
    @test MetabolicPathway((:A_ext,), ((:Enz1, (:A_ext,), (:B,)),)) isa MetabolicPathway
    @test MetabolicPathway{(:A_ext,),((:Enz1, (:A_ext,), (:B,)),)}() isa MetabolicPathway
    @test MetabolicPathway((:A_ext,), ((:Enz1, (:A_ext,), (:B,)),)) ===
          MetabolicPathway{(:A_ext,),((:Enz1, (:A_ext,), (:B,)),)}()
    @test_throws ErrorException MetabolicPathway(:A_ext, ((:Enz1, (:A_ext,), (:B,)),))
    @test_throws ErrorException MetabolicPathway((:A_ext,), (:Enz1, (:A_ext,), (:B,)))

    # Test MetabolicPathway functions
    test_pathway = MetabolicPathway(
        (:A_media,),
        (
            (:Enz1, (:A_media,), (:A,), (:Activator1,), (:Inhibitor1,)),
            (:Enz2, (:A,), (:B, :B), (), (:Inhibitor2, :Inhibitor3)),
            (:Enz3, (:B,), (:C,), (:Activator2, :Activator3), ()),
            (:Enz4, (:C, :C), (:D,), (:Activator4,), (:Inhibitor4,)),
        ),
    )

    test_pathway_metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0, D = 1.0)
    expected_constant_metabolites = (:A_media,)
    expected_enzyme_names = (:Enz1, :Enz2, :Enz3, :Enz4)
    expected_substrate_names = ((:A_media,), (:A,), (:B,), (:C, :C))
    expected_product_names = ((:A,), (:B, :B), (:C,), (:D,))
    expected_activator_names =
        ((:Activator1,), (), (:Activator2, :Activator3), (:Activator4,))
    expected_inhibitor_names =
        ((:Inhibitor1,), (:Inhibitor2, :Inhibitor3), (), (:Inhibitor4,))
    expected_all_metabolite_names = (
        :A_media,
        :A,
        :B,
        :C,
        :D,
        :Activator1,
        :Activator2,
        :Activator3,
        :Activator4,
        :Inhibitor1,
        :Inhibitor2,
        :Inhibitor3,
        :Inhibitor4,
    )
    expected_stoichiometric_matrix = [
        0 0 0 0
        1 -1 0 0
        0 2 -1 0
        0 0 1 -2
        0 0 0 1
    ]
    test_params = LVector(Enz1_Keq = 10.0, Enz2_Keq = 5.0, Enz3_Keq = 2.0, Enz4_Keq = 1.5)

    # Test metabolite_names function
    @test constant_metabolites(test_pathway) == expected_constant_metabolites
    constant_metabolites(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark constant_metabolites($test_pathway)
    @test minimum(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test enzymes function
    @test enzymes(test_pathway) == expected_enzyme_names
    enzymes(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark enzymes($test_pathway)
    # Non-@generated accessor: does real runtime work (_generate_Enzymes + map),
    # so wall-time tracks the CI CPU (~17 ns on AMD znver3, ~42 ns on Intel icelake).
    # 50 ns gives margin above the slowest observed runner. allocs == 0 is the exact guard.
    @test minimum(benchmark_result.times) <= 50 #ns
    @test benchmark_result.allocs == 0

    # Test substrates function
    @test substrates(test_pathway) == expected_substrate_names
    substrates(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark substrates($test_pathway)
    @test minimum(benchmark_result.times) <= 50 #ns (non-@generated; see enzymes note)
    @test benchmark_result.allocs == 0

    # Test products function
    @test products(test_pathway) == expected_product_names
    products(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark products($test_pathway)
    @test minimum(benchmark_result.times) <= 50 #ns (non-@generated; see enzymes note)
    @test benchmark_result.allocs == 0

    # Test stoichiometric_matrix function
    @test stoichiometric_matrix(test_pathway) == expected_stoichiometric_matrix
    stoichiometric_matrix(test_pathway) isa Matrix{Int}
    benchmark_result = @benchmark stoichiometric_matrix($test_pathway)
    @test minimum(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test activators function
    @test activators(test_pathway) == expected_activator_names
    activators(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark activators($test_pathway)
    @test minimum(benchmark_result.times) <= 50 #ns (non-@generated; see enzymes note)
    @test benchmark_result.allocs == 0

    # Test inhibitors function
    @test inhibitors(test_pathway) == expected_inhibitor_names
    inhibitors(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark inhibitors($test_pathway)
    @test minimum(benchmark_result.times) <= 50 #ns (non-@generated; see enzymes note)
    @test benchmark_result.allocs == 0

    enzyme_instances = CellMetabolismBase._generate_Enzymes(test_pathway)
    @test name(enzyme_instances[1]) == :Enz1
    @test name(enzyme_instances[2]) == :Enz2
    @test substrates(enzyme_instances[1]) == (:A_media,)
    @test substrates(enzyme_instances[4]) == (:C, :C)
    @test products(enzyme_instances[1]) == (:A,)
    @test products(enzyme_instances[2]) == (:B, :B)
    @test activators(enzyme_instances[1]) == (:Activator1,)
    @test activators(enzyme_instances[2]) == ()
    @test inhibitors(enzyme_instances[1]) == (:Inhibitor1,)
    @test inhibitors(enzyme_instances[3]) == ()
    @test enzymes(test_pathway) == map(name, enzyme_instances)
    @test substrates(test_pathway) == map(substrates, enzyme_instances)
    @test products(test_pathway) == map(products, enzyme_instances)
    @test activators(test_pathway) == map(activators, enzyme_instances)
    @test inhibitors(test_pathway) == map(inhibitors, enzyme_instances)

    single_ratios = map(
        enzyme -> disequilibrium_ratio(enzyme, test_pathway_metabs, test_params),
        enzyme_instances,
    )
    @test disequilibrium_ratios(test_pathway, test_pathway_metabs, test_params) ==
          single_ratios
    @test disequilibrium_ratio(enzyme_instances[1], test_pathway_metabs, test_params) ≈
          (test_pathway_metabs.A / test_pathway_metabs.A_media) / test_params.Enz1_Keq

    params_missing = LVector(Enz1_Keq = 10.0, Enz2_Keq = 5.0, Enz4_Keq = 1.5)
    @test_throws ArgumentError disequilibrium_ratio(
        enzyme_instances[3],
        test_pathway_metabs,
        params_missing,
    )

    # Test reactants function
    @test reactants(test_pathway) == (:A_media, :A, :B, :C, :D)
    reactants(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark reactants($test_pathway)
    @test minimum(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test metabolites function
    @test metabolites(test_pathway) == expected_all_metabolite_names
    metabolites(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark metabolites($test_pathway)
    @test minimum(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    shorthand_pathway = MetabolicPathway(
        (:A_media,),
        ((:Simple, (:A_media,), (:A,)), (:Regulated, (:A,), (:B,), (), ())),
    )
    @test activators(shorthand_pathway) == ((), ())
    @test inhibitors(shorthand_pathway) == ((), ())
    shorthand_enzymes = CellMetabolismBase._generate_Enzymes(shorthand_pathway)
    @test activators(shorthand_enzymes[1]) == ()
    @test inhibitors(shorthand_enzymes[1]) == ()
    @test activators(shorthand_pathway) == map(activators, shorthand_enzymes)
    @test inhibitors(shorthand_pathway) == map(inhibitors, shorthand_enzymes)

    @test_throws ErrorException rate(Enzyme(:Random, (:X,), (:Y,)), 1.0, 2.0)
end

@testitem "remove_regulation throws error for undefined enzyme" begin
    using CellMetabolismBase
    using LabelledArrays

    enzyme = CellMetabolismBase.Enzyme(:Demo, (:S,), (:P,), (:Act1, :Act2), (:Inh1,))
    params = LVector(
        Demo_K_a_Act1 = 1.0,
        Demo_K_a_Act2 = 1.5,
        Demo_K_i_Inh1 = 2.0,
        Demo_L = 3.0,
        Demo_Vmax = 5.0,
    )

    @test_throws ErrorException CellMetabolismBase.remove_regulation(enzyme, params)
    @test_throws ErrorException CellMetabolismBase.remove_regulation(enzyme, params, Val(:Act1))
end

@testitem "remove_regulation user extension example" begin
    using CellMetabolismBase
    using LabelledArrays

    enzyme = CellMetabolismBase.Enzyme(:Custom, (:S,), (:P,), (:Act,), ())
    params = LVector(
        Custom_K_a_Act = 2.0,
        Custom_L = 4.0,
        Custom_bias = 10.0,
    )

    # User extends the function for specific regulator
    function CellMetabolismBase.remove_regulation(
        enzyme::CellMetabolismBase.Enzyme{:Custom,(:S,),(:P,),(:Act,),()},
        params,
        ::Val{:Act},
    )
        params = deepcopy(params)
        setproperty!(params, :Custom_K_a_Act, Inf)
        params.Custom_bias = 0.0
        return params
    end

    # User extends the function for removing all regulation
    function CellMetabolismBase.remove_regulation(
        enzyme::CellMetabolismBase.Enzyme{:Custom,(:S,),(:P,),(:Act,),()},
        params,
    )
        params = deepcopy(params)
        setproperty!(params, :Custom_L, 0.0)
        params = CellMetabolismBase.remove_regulation(enzyme, params, Val(:Act))
        return params
    end

    # Test single regulator removal
    act_removed = CellMetabolismBase.remove_regulation(enzyme, params, Val(:Act))
    @test isinf(act_removed.Custom_K_a_Act)
    @test act_removed.Custom_L == params.Custom_L
    @test act_removed.Custom_bias == 0.0

    # Test all regulation removal
    neutral_params = CellMetabolismBase.remove_regulation(enzyme, params)
    @test isinf(neutral_params.Custom_K_a_Act)
    @test neutral_params.Custom_L == 0.0
    @test neutral_params.Custom_bias == 0.0
end

@testsnippet GenerateRandomStoichiometricMatrix begin
    using SparseArrays, Distributions, StatsBase, Random

    """
    Helper to calculate the probability weights that determine the number of reactions a
    metabolite is involved in according to a power law.
    """
    function get_powerlaw_weights(m::Int, alpha::Float64)
        return [1.0 / (i^alpha) for i = 1:m]
    end

    """
    Helper to generate a random stoichiometric matrix that has a scale-free network topology
    (Jeon et. al. 2000). Each reaction has a minimum of 2 participants and an average of
    `avg_participants`. 90% of the coefficients in the network are -1 or 1, 10% of the
    coefficients are -2 or 2. The steepness of the power law that defines the scale-free
    network is set by `alpha` (default value matches finding of Jeon et. al. 2000).

    References
        Jeong, Hawoong, et al. "The large-scale organization of metabolic networks." Nature
        407.6804 (2000): 651-654.
    """
    function generate_random_stoichiometric_matrix(
        m::Int,
        n::Int;
        avg_participants::Float64 = 4.0,
        alpha::Float64 = 0.833334,
    )
        weights = get_powerlaw_weights(m, alpha)

        # Construct vector of weight values
        prob_weights = Weights(weights)

        # Shuffle the indices so the mega-hubs aren't always rows 1, 2, and 3
        metabolite_ids = randperm(m)

        I = Int[] # row index
        J = Int[] # column index
        V = Int[] # value

        rxn_size_dist = Poisson(avg_participants - 2.0)

        for j = 1:n
            num_participants = min(m, 2.0 + rand(rxn_size_dist))

            participants = Set{Int}()
            while length(participants) < num_participants
                # Draw a metabolite strictly according to the Power Law
                chosen_index = sample(metabolite_ids, prob_weights)
                push!(participants, chosen_index)
            end

            for i in participants
                coeff = rand() < 0.9 ? rand([-1, 1]) : rand([-2, 2])
                push!(I, i), push!(J, j), push!(V, coeff)
            end
        end

        return sparse(I, J, V, m, n)
    end
end

@testitem "Random stoichiometric matrix generation" setup =
    [GenerateRandomStoichiometricMatrix] begin
    using Random

    # Create a random matrix for all following tests
    m = rand(10:50)
    n = rand(10:50)
    random_S = generate_random_stoichiometric_matrix(m, n)

    @testset "The dimensions and type are correct" begin
        @test size(random_S, 1) == m && size(random_S, 2) == n
        @test random_S isa SparseMatrixCSC
    end

    @testset "The desity of the matrix is determined by the `avg_participants`" begin
        avg_participants = 4.0
        expected_density = avg_participants / m
        actual_density = nnz(random_S) / (m * n)
        @test isapprox(expected_density, actual_density, rtol = 0.1)
    end

    @testset "Coefficient values must only be -2, -1, 1, or 2" begin
        valid_coeffs = Set([-2, -1, 1, 2])
        @test all(v -> v in valid_coeffs, nonzeros(random_S))
    end

    @testset "Reaction sizes should be ≥ 2 or ≤ `m`" begin
        # Get the number of non-zeros in each column
        col_participant_counts = diff(random_S.colptr)
        @test minimum(col_participant_counts) >= 2
        @test maximum(col_participant_counts) <= m
    end

    @testset "Scale-free network topology" begin
        # In a scale-free network, mega-hub metabolites should be involved in many more
        # reactions than the median. In a normal/Poisson network, the metabolite involved
        # in the maximum number of reactions is usually close to the median
        Random.seed!(42)
        S = generate_random_stoichiometric_matrix(1000, 5000)

        # Get degrees of all metabolites
        row_degrees = zeros(Int, 1000)
        for r in rowvals(S)
            row_degrees[r] += 1
        end

        max_degree = maximum(row_degrees)
        median_degree = median(row_degrees)

        @test max_degree > (median_degree * 10)
    end

    @testset "Power Law Weight Generation" begin
        m = 5
        alpha = 0.8333

        # Generate the weights
        weights = get_powerlaw_weights(m, alpha)

        # 1. The first weight should always be exactly 1.0 (1 / 1^alpha)
        @test weights[1] == 1.0

        # 2. Check the mathematical decline of the weights
        # For element 2, weight should be 1 / (2^0.8333) ≈ 0.56123
        @test isapprox(weights[2], 1.0 / (2^alpha), atol = 1e-5)

        # For element 5, weight should be 1 / (5^0.8333) ≈ 0.26116
        @test isapprox(weights[5], 1.0 / (5^alpha), atol = 1e-5)

        # 3. Ensure all weights are strictly positive and descending
        @test all(weights .> 0)
        @test issorted(weights, rev = true)
    end
end

@testitem "Conserved Moieties Tests" setup = [GenerateRandomStoichiometricMatrix] begin
    @testset "Fourier-Motzkin Fundamentals" begin
        # Test 10 random stoichiometric matrices to ensure fundamental checks pass
        i = 1
        number_of_tests = 10

        # Store the results of each test
        nullspace_test_results = Bool[]
        positive_test_results = Bool[]
        metabolites_present_result = Bool[]
        non_zero_moieties_result = Bool[]
        simplest_form_result = Bool[]
        minimal_pools_result = Bool[]

        keep_testing = true
        while keep_testing
            reactions = rand(1:100)
            # 1.3 times more metabolites than reactions (similar to Glycolysis and TCA cycle)
            metabolites = round(Int, 1.31 * reactions)
            random_matrix = generate_random_stoichiometric_matrix(metabolites, reactions)
            R_random = CellMetabolismBase._fourier_motzkin(random_matrix)

            # Test 1: The solution is the left null space of the stoichiometric matrix
            push!(nullspace_test_results, iszero(R_random * random_matrix))

            # Test 2: There are no negative values in the matrix of conserved moieties
            push!(positive_test_results, all(R_random .>= 0))

            # Test 3: All metabolites in the stoichiometric matrix are present in the
            # conserved moieties matrix
            push!(metabolites_present_result, size(R_random, 2) == metabolites)

            for row_i in eachrow(R_random)

                # Test 4: No conserved moiety is entirely zero
                push!(non_zero_moieties_result, !iszero(row_i))

                # Test 5: Each conserved moiety should be in its simplest form, so the
                # greatest common divisor should be equal to 1
                push!(simplest_form_result, gcd(row_i) == 1)
                for row_j in eachrow(R_random)
                    if row_i != row_j

                        # Test 6. None of the conserved moieties should be a subset of
                        # another conserved moiety
                        push!(
                            minimal_pools_result,
                            !issubset(
                                findall(x -> x != 0, row_i),
                                findall(x -> x != 0, row_j),
                            ),
                        )
                    end
                end
            end
            if i == number_of_tests
                keep_testing = false
            end
            i += 1
        end

        @test all(nullspace_test_results)
        @test all(positive_test_results)
        @test all(metabolites_present_result)
        @test all(non_zero_moieties_result)
        @test all(simplest_form_result)
        @test all(minimal_pools_result)
    end

    @testset "_find_best_column returns the best reaction to process next" begin
        # Reactions
        # 1: ATP + B + C ⇾ D + E + ADP
        # 2: ATP + D + E ⇾ F + G + ADP + H
        # 3: I ⇾ J
        # 4: K ⇾ L

        # Costs
        # Cost formula: (num_pos * num_neg) - (num_pos + num_neg)
        # 1: 3*3 - (3+3) = 9 - 6 = 3
        # 2: 3*4 - (3+4) = 12 - 7 = 5
        # 3: 1*1 - (1-1) = 0
        # 4: 1*1 - (1-1) = 0
        S = [
            -1 -1 0 0 # ATP
            -1 0 0 0 # B
            -1 0 0 0 # C
            1 -1 0 0 # D
            1 -1 0 0 # E
            0 1 0 0 # F
            0 1 0 0 # G
            1 1 0 0 # ADP
            0 1 0 0 # H
            0 0 1 0 # I
            0 0 -1 0 # J
            0 0 0 1 # K
            0 0 0 -1 # L
        ]

        @testset "Combinatorial explosion prevention" begin
            # Since ATP and ADP are involved in reactions 1 and 2, the best reaction to
            # process first is 3 because it will cause the least combinatorial explosion of
            # terms and it appears before reaction 4 (which has the same cost)
            best_idx = CellMetabolismBase._find_best_column(S, 1:4)
            @test best_idx == 3
        end

        @testset "Tie breaker behaviour" begin
            # If evaluated on [3, 4], both have cost 0.
            # The algorithm strictly uses `< min_cost`, so it returns the first evaluated.
            @test CellMetabolismBase._find_best_column(S, [3, 4]) == 3
            @test CellMetabolismBase._find_best_column(S, [4, 3]) == 4
        end

        @testset "Subset of remaining columns" begin
            # If column 3 is removed, it should evaluate [1, 2, 4]
            # Costs are: C1 => 3, C2 => 5, C4 => 0. Best is C4.
            @test CellMetabolismBase._find_best_column(S, [1, 2, 4]) == 4

            # If column 3 and 4 are removed, it should evaluate [1, 2]
            # Costs are: C1 => 3, C2 => 5. Best is C1.
            @test CellMetabolismBase._find_best_column(S, [1, 2]) == 1
        end

        @testset "Empty remaining columns" begin
            # The function initializes `best_col` to -1.
            # If passed an empty list, it should return -1 safely without bounds errors.
            @test CellMetabolismBase._find_best_column(S, Int[]) == -1
        end
    end

    @testset "Toy networks tests" begin
        @testset "Simple conversion: 2A -> 2B" begin
            S1 = [
                -2; 2;;
            ]
            R1 = CellMetabolismBase._fourier_motzkin(S1)
            @test size(R1, 1) == 1 # 1 conserved pool
            @test size(R1, 2) == 2 # 2 metabolites
            # The pool should be 1*A + 1*B
            @test R1[1, 1] == R1[1, 2] && R1[1, :] == [1, 1]
            @test iszero(R1 * S1)
            @test all(R1 .>= 0)
        end

        @testset "Linear Chain: A -> B -> C" begin
            S2 = [
                -1 0
                1 -1
                0 1
            ]
            R2 = CellMetabolismBase._fourier_motzkin(S2)
            @test size(R2, 1) == 1 # 1 conserved pool
            @test size(R2, 2) == 3 # 3 metabolites
            # The pool should be 1*A + 1*B + 1*C
            @test R2[1, 1] == R2[1, 2] == R2[1, 3] && R2[1, :] == [1, 1, 1]
            @test iszero(R2 * S2)
            @test all(R2 .>= 0)
        end

        @testset "Open System: -> A ->" begin
            S3 = [
                1 -1
            ]
            R3 = CellMetabolismBase._fourier_motzkin(S3)
            @test size(R3, 1) == 0 # 0 conserved pools (nothing is conserved)
            @test size(R3, 2) == 1 # 1 metabolite
            @test iszero(R3 * S3)
        end

        @testset "Disconnected metabolites" begin
            S4 = [
                0; 0;;
            ]
            R4 = CellMetabolismBase._fourier_motzkin(S4)
            # 2 conserved pools (each metabolite independently)
            @test size(R4, 1) == 2 # 2 conserved pools
            @test size(R4, 2) == 2 # 2 metabolites
            @test R4[1, :] == [1, 0] && R4[2, :] == [0, 1]
            @test iszero(R4 * S4)
            @test all(R4 .>= 0)
        end

        @testset "Closed Loop: A → B → C → A" begin
            S5 = [
                -1 0 1
                1 -1 0
                0 1 -1
            ]
            R5 = CellMetabolismBase._fourier_motzkin(S5)
            # One conserved moiety (1*A + 1*B + 1*C) due to pathway using and making A
            @test size(R5, 1) == 1
            @test size(R5, 2) == 3 # 3 metabolites
            @test R5[1, :] == [1, 1, 1]
            @test iszero(R5 * S5)
            @test all(R5 .>= 0)
        end

        @testset "Carrier Cofactor: A + ATP → B + ADP" begin
            S6 = [
                -1; -1; 1; 1;;
            ]
            R6 = CellMetabolismBase._fourier_motzkin(S6)
            # Four conserved pools: 1*A + 1*B, 1*B + 1*ATP, 1*A + 1*ADP, 1*B + 1*ATP
            @test size(R6, 1) == 4
            @test size(R6, 2) == 4 # 4 metabolites
            # 1*A + 1*B
            @test R6[1, 1] == R6[1, 3] && R6[1, :] == [1, 0, 1, 0]
            # 1*B + 1*ATP
            @test R6[2, 2] == R6[2, 3] && R6[2, :] == [0, 1, 1, 0]
            # 1*A + 1*ADP
            @test R6[3, 1] == R6[3, 4] && R6[3, :] == [1, 0, 0, 1]
            # 1*ATP + 1*ADP
            @test R6[4, 2] == R6[4, 4] && R6[4, :] == [0, 1, 0, 1]
            @test iszero(R6 * S6)
            @test all(R6 .>= 0)
        end

        @testset "Branched Pathway: A → B, A → C" begin
            S7 = [
                -1 -1
                1 0
                0 1
            ]
            R7 = CellMetabolismBase._fourier_motzkin(S7)
            # One conserved pool (1*A + 1*B + 1*C)
            @test size(R7, 1) == 1
            @test size(R7, 2) == 3 # 3 metabolites
            @test R7[1, 1] == R7[1, 2] == R7[1, 3] && R7[1, :] == [1, 1, 1]
            @test iszero(R7 * S7)
            @test all(R7 .>= 0)
        end

        @testset "Disconnected pathways: A → B → C, D → E → F" begin
            S8 = [
                -1 0 0 0
                1 -1 0 0
                0 1 0 0
                0 0 -1 0
                0 0 1 -1
                0 0 0 1
            ]
            R8 = CellMetabolismBase._fourier_motzkin(S8)

            # Two conserved pools: 1*A + 1*B + 1*C and 1*D + 1*E + 1*F
            @test size(R8, 1) == 2
            @test size(R8, 2) == 6 # 6 metabolites
            @test R8[1, 1] == R8[1, 2] == R8[1, 3] && R8[1, :] == [1, 1, 1, 0, 0, 0]
            @test R8[2, 4] == R8[2, 5] == R8[2, 6] && R8[2, :] == [0, 0, 0, 1, 1, 1]

            @test iszero(R8 * S8)
            @test all(R8 .>= 0)
        end

        @testset "Different stoichiometries: 1A + 2B ⇾ 2C" begin
            S9 = [
                -1; -2; 2;;
            ]
            R9 = CellMetabolismBase._fourier_motzkin(S9)
            # Two conserved pools: 2*A + 1*C and 1*B + 1*C
            @test size(R9, 1) == 2
            @test size(R9, 2) == 3 # 3 metabolites
            @test R9[1, 1] == R9[1, 3] * 2 && R9[2, 2] == R9[2, 3]
            @test R9[1, :] == [2, 0, 1] && R9[2, :] == [0, 1, 1]
            @test iszero(R9 * S9)
            @test all(R9 .>= 0)
        end
    end

    @testset "conserved_moieties is correct for a known pathway" begin
        glycolysis_pathway = MetabolicPathway(
            (:Glucose_media, :Lactate_media),
            (
                (:GLUT, (:Glucose_media,), (:Glucose,)),
                (:HK1, (:Glucose, :ATP), (:G6P, :ADP)),
                (:GPI, (:G6P,), (:F6P,)),
                (:PFKP, (:F6P, :ATP), (:F16BP, :ADP)),
                (:ALDO, (:F16BP,), (:GAP, :DHAP)),
                (:TPI, (:DHAP,), (:GAP,)),
                (:GAPDH, (:GAP, :Pi, :NAD), (:BPG, :NADH)),
                (:PGK, (:BPG, :ADP), (:Three_PG, :ATP)),
                (:PGM, (:Three_PG,), (:Two_PG,)),
                (:ENO, (:Two_PG,), (:PEP,)),
                (:PKM, (:PEP, :ADP), (:Pyruvate, :ATP)),
                (:LDH, (:Pyruvate, :NADH), (:Lactate, :NAD)),
                (:MCT, (:Lactate,), (:Lactate_media,)),
                (:ATPase, (:ATP,), (:ADP, :Pi)),
                (:AK, (:ADP, :ADP), (:ATP, :AMP)),
            ),
        )
        expected_glycolysis_conserved_moieties = join(
            [
                "1: 1⋅Glucose_media",
                "2: 1⋅Lactate_media",
                "3: 1⋅NAD + 1⋅NADH",
                "4: 1⋅ATP + 1⋅ADP + 1⋅AMP",
                "5: 1⋅NAD + 1⋅BPG + 1⋅Three_PG + 1⋅Two_PG + 1⋅PEP + 1⋅Pyruvate",
                "6: 2⋅ATP + 1⋅G6P + 1⋅ADP + 1⋅F6P + 2⋅F16BP + 1⋅GAP + 1⋅DHAP + 1⋅Pi + " *
                "2⋅BPG + 1⋅Three_PG + 1⋅Two_PG + 1⋅PEP",
            ],
            "\n",
        )
        glycolysis_conserved_moieties = conserved_moieties(glycolysis_pathway)
        @test glycolysis_conserved_moieties == expected_glycolysis_conserved_moieties
    end
end
