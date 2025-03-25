@testitem "Accessor functions for MetabolicPathways" begin
    using BenchmarkTools, LabelledArrays

    # Test basic Enzyme constructor and its errors
    @test Enzyme(:Enz1, (:A,), (:B,)) isa Enzyme
    @test Enzyme(:Enz1,(:A,),(:B,)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,)) === Enzyme{:Enz1,(:A,),(:B,),(),()}()
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1,), (:Inh1,)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1, :Act2), (:Inh1, :Inh2)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (), (:Inh1,)) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1,), ()) isa Enzyme
    @test Enzyme(:Enz1, (:A,), (:B,), (:Act1,), (:Inh1,)) === Enzyme{:Enz1,(:A,),(:B,),(:Act1,),(:Inh1,)}()
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
    expected_activator_names = ((:Activator1,), (), (:Activator2, :Activator3), (:Activator4,))
    expected_inhibitor_names = ((:Inhibitor1,), (:Inhibitor2, :Inhibitor3), (), (:Inhibitor4,))
    expected_stoichiometric_matrix = [
        0 0 0 0
        1 -1 0 0
        0 2 -1 0
        0 0 1 -2
        0 0 0 1
    ]

    @test constant_metabs(test_pathway) == expected_const_metabs
    constant_metabs(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark constant_metabs($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    @test enzyme_names(test_pathway) == expected_enzyme_names
    enzyme_names(test_pathway) isa Tuple{Vararg{Symbol}}
    benchmark_result = @benchmark enzyme_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    @test substrate_names(test_pathway) == expected_substrate_names
    substrate_names(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark substrate_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    @test product_names(test_pathway) == expected_product_names
    product_names(test_pathway) isa Tuple{Vararg{Tuple{Vararg{Symbol}}}}
    benchmark_result = @benchmark product_names($test_pathway)
    @test mean(benchmark_result.times) <= 10 #ns
    @test benchmark_result.allocs == 0

    @test stoichiometric_matrix(test_pathway, test_pathway_metabs) ==
          expected_stoichiometric_matrix
    stoichiometric_matrix(test_pathway, test_pathway_metabs) isa Matrix{Int}
    benchmark_result = @benchmark stoichiometric_matrix($test_pathway, $test_pathway_metabs)
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

    @test_throws ErrorException enzyme_rate(Enzyme(:Random,(:X,),(:Y,)), 1.0, 2.0)
end
