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
    @test mean(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test enzymes function
    @test enzymes(test_pathway) == expected_enzyme_names
    enzymes(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark enzymes($test_pathway)
    @test mean(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test substrates function
    @test substrates(test_pathway) == expected_substrate_names
    substrates(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark substrates($test_pathway)
    @test mean(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test products function
    @test products(test_pathway) == expected_product_names
    products(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark products($test_pathway)
    @test mean(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test stoichiometric_matrix function
    @test stoichiometric_matrix(test_pathway) == expected_stoichiometric_matrix
    stoichiometric_matrix(test_pathway) isa Matrix{Int}
    benchmark_result = @benchmark stoichiometric_matrix($test_pathway)
    @test mean(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test activators function
    @test activators(test_pathway) == expected_activator_names
    activators(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark activators($test_pathway)
    @test mean(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test inhibitors function
    @test inhibitors(test_pathway) == expected_inhibitor_names
    inhibitors(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark inhibitors($test_pathway)
    @test mean(benchmark_result.times) <= 20 #ns
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
    @test mean(benchmark_result.times) <= 20 #ns
    @test benchmark_result.allocs == 0

    # Test metabolites function
    @test metabolites(test_pathway) == expected_all_metabolite_names
    metabolites(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark metabolites($test_pathway)
    @test mean(benchmark_result.times) <= 20 #ns
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
