@testitem "Validation of MetabolicPathway" begin
    using LabelledArrays, BenchmarkTools, OrdinaryDiffEq

    test_pathway = MetabolicPathway(
        (:A_media,),
        (
            (:Enz1, (:A_media,), (:A,)),
            (:Enz2, (:A,), (:B, :B)),
            (:Enz3, (:B,), (:C,)),
            (:Enz4, (:C, :C), (:D,)),
        ),
    )

    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:Enz1,(:A_media,),(:A,)},
        metabs,
        params,
    )
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

    metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0, D = 1.0)
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


    CellMetabolismBase.validate_enzymes(test_pathway, metabs, params)

    CellMetabolismBase.validate_metabolic_pathway(
        test_pathway,
        metabs,
        params,
    )
end
