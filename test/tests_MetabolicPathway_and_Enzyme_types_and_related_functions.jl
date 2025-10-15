@testitem "Accessor functions for MetabolicPathways" begin
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
    expected_const_metabs = (:A_media,)
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
    test_params = LVector(
        Enz1_Keq = 10.0,
        Enz2_Keq = 5.0,
        Enz3_Keq = 2.0,
        Enz4_Keq = 1.5,
    )

    # Test metabolite_names function
    @test constant_metabs(test_pathway) == expected_const_metabs
    constant_metabs(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark constant_metabs($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    # Test enzyme_names function
    @test enzyme_names(test_pathway) == expected_enzyme_names
    enzyme_names(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark enzyme_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    # Test substrate_names function
    @test substrate_names(test_pathway) == expected_substrate_names
    substrate_names(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark substrate_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    # Test product_names function
    @test product_names(test_pathway) == expected_product_names
    product_names(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark product_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    # Test stoichiometric_matrix function
    @test stoichiometric_matrix(test_pathway) == expected_stoichiometric_matrix
    stoichiometric_matrix(test_pathway) isa Matrix{Int}
    benchmark_result = @benchmark stoichiometric_matrix($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    # Test activator_names function
    @test activator_names(test_pathway) == expected_activator_names
    activator_names(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark activator_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    # Test inhibitor_names function
    @test inhibitor_names(test_pathway) == expected_inhibitor_names
    inhibitor_names(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark inhibitor_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    enzymes = CellMetabolismBase._generate_Enzymes(test_pathway)
    @test enzyme_name(enzymes[1]) == :Enz1
    @test enzyme_name(enzymes[2]) == :Enz2
    @test substrates_name(enzymes[1]) == (:A_media,)
    @test substrates_name(enzymes[4]) == (:C, :C)
    @test products_name(enzymes[1]) == (:A,)
    @test products_name(enzymes[2]) == (:B, :B)
    @test activators_name(enzymes[1]) == (:Activator1,)
    @test activators_name(enzymes[2]) == ()
    @test inhibitors_name(enzymes[1]) == (:Inhibitor1,)
    @test inhibitors_name(enzymes[3]) == ()
    @test enzyme_names(test_pathway) == map(enzyme_name, enzymes)
    @test substrate_names(test_pathway) == map(substrates_name, enzymes)
    @test product_names(test_pathway) == map(products_name, enzymes)
    @test activator_names(test_pathway) == map(activators_name, enzymes)
    @test inhibitor_names(test_pathway) == map(inhibitors_name, enzymes)

    single_ratios = map(enzyme -> disequilibrium_ratio(enzyme, test_pathway_metabs, test_params), enzymes)
    @test disequilibrium_ratios(test_pathway, test_pathway_metabs, test_params) == single_ratios
    @test disequilibrium_ratio(enzymes[1], test_pathway_metabs, test_params) ≈
          (test_pathway_metabs.A / test_pathway_metabs.A_media) / test_params.Enz1_Keq

    params_missing = LVector(
        Enz1_Keq = 10.0,
        Enz2_Keq = 5.0,
        Enz4_Keq = 1.5,
    )
    @test_throws ArgumentError disequilibrium_ratio(enzymes[3], test_pathway_metabs, params_missing)

    # Test reactant_names function
    @test reactant_names(test_pathway) == (:A_media, :A, :B, :C, :D)
    reactant_names(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark reactant_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    # Test all_metabolite_names function
    @test all_metabolite_names(test_pathway) == expected_all_metabolite_names
    all_metabolite_names(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark all_metabolite_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    @test_throws ErrorException enzyme_rate(Enzyme(:Random, (:X,), (:Y,)), 1.0, 2.0)
end
