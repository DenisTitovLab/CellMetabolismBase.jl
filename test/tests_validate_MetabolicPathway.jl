@testitem "MetabolicPathway Validation Tests" setup=[TestMetabolicPathway] begin
    using LabelledArrays, BenchmarkTools, OrdinaryDiffEqFIRK

    # Define a common test pathway with different enzyme types
    test_pathway = TestMetabolicPathway.test_pathway

    # Define standard metabolites and parameters
    valid_metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0, D = 1.0)
    valid_params = LVector(
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

    # ------- Test: Full validation with standard pathway -------
    # Test that the overall validation function works
    @test_nowarn CellMetabolismBase.validate_MetabolicPathway(
        test_pathway,
        valid_metabs,
        valid_params,
    )
    # Test the enzyme validation function works
    @test_nowarn CellMetabolismBase.validate_enzyme_rates(
        test_pathway,
        valid_metabs,
        valid_params,
    )

    # ------- Test: Check metabolites found in initial conditions -------
    # Test 1: Missing metabolite - should fail
    missing_metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0) # Missing C and D
    @test_throws ErrorException CellMetabolismBase.validate_MetabolicPathway(
        test_pathway,
        missing_metabs,
        valid_params,
    )

    # Test 2: Extra metabolite - should pass
    extra_metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0, D = 1.0, E = 0.5) # Extra E
    @test_nowarn CellMetabolismBase.validate_MetabolicPathway(
        test_pathway,
        extra_metabs,
        valid_params,
    )

    # ------- Test: enzyme_rate() is defined and returns Real number -------
    # Test with a simple pathway
    simple_pathway = MetabolicPathway((:A_ext,), ((:SimpleEnz, (:A_ext,), (:B,)),))
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:SimpleEnz,(:A_ext,),(:B,)},
        metabs,
        params,
    )
        return params.SimpleEnz_Vmax * (metabs.A_ext - metabs.B / params.SimpleEnz_Keq) / (
            params.SimpleEnz_Km *
            (1 + metabs.A_ext / params.SimpleEnz_Km + metabs.B / params.SimpleEnz_Ki)
        )
    end
    simple_metabs = LVector(A_ext = 1.0, B = 0.5)
    simple_params = LVector(
        SimpleEnz_Vmax = 1.0,
        SimpleEnz_Km = 0.1,
        SimpleEnz_Ki = 0.5,
        SimpleEnz_Keq = 2.0,
    )
    # Should not throw an error
    @test_nowarn CellMetabolismBase.validate_enzyme_rates(
        simple_pathway,
        simple_metabs,
        simple_params,
    )

    # Test with a non-Real return type (String)
    simple_string_pathway =
        MetabolicPathway((:A_ext,), ((:SimpleEnzString, (:A_ext,), (:B,)),))
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:SimpleEnzString,(:A_ext,),(:B,)},
        metabs,
        params,
    )
        return "Not a Real number"
    end
    simple_string_params = LVector(
        SimpleEnzString_Vmax = 1.0,
        SimpleEnzString_Km = 0.1,
        SimpleEnzString_Ki = 0.5,
        SimpleEnzString_Keq = 2.0,
    )
    @test_throws ErrorException CellMetabolismBase.validate_enzyme_rates(
        simple_string_pathway,
        simple_metabs,
        simple_string_params,
    )

    # Test with a non-Real return type (Array)
    simple_array_pathway =
        MetabolicPathway((:A_ext,), ((:SimpleEnzArray, (:A_ext,), (:B,)),))
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:SimpleEnzArray,(:A_ext,),(:B,)},
        metabs,
        params,
    )
        return [1.0, 2.0]  # Array instead of Real
    end
    simple_array_params = LVector(
        SimpleEnzArray_Vmax = 1.0,
        SimpleEnzArray_Km = 0.1,
        SimpleEnzArray_Ki = 0.5,
        SimpleEnzArray_Keq = 2.0,
    )
    @test_throws ErrorException CellMetabolismBase.validate_enzyme_rates(
        simple_array_pathway,
        simple_metabs,
        simple_array_params,
    )

    # Test with a non-Real return type (Nothing)
    simple_nothing_pathway =
        MetabolicPathway((:A_ext,), ((:SimpleEnzNothing, (:A_ext,), (:B,)),))
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:SimpleEnzNothing,(:A_ext,),(:B,)},
        metabs,
        params,
    )
        return nothing  # Nothing instead of Real
    end
    simple_nothing_params = LVector(
        SimpleEnzNothing_Vmax = 1.0,
        SimpleEnzNothing_Km = 0.1,
        SimpleEnzNothing_Ki = 0.5,
        SimpleEnzNothing_Keq = 2.0,
    )
    @test_throws ErrorException CellMetabolismBase.validate_enzyme_rates(
        simple_nothing_pathway,
        simple_metabs,
        simple_nothing_params,
    )

    # Missing enzyme rate function test
    missing_function_pathway = MetabolicPathway((:X,), ((:MissingFunc, (:X,), (:Y,)),))

    missing_metabs = LVector(X = 1.0, Y = 0.5)
    missing_params = LVector(MissingFunc_Vmax = 1.0)

    # Should throw an error when enzyme_rate function is not defined
    @test_throws Exception CellMetabolismBase.validate_enzyme_rates(
        missing_function_pathway,
        missing_metabs,
        missing_params,
    )


    # ------- Test: Enzyme rate is positive when products absent, substrates present -------
    # Create a specialized test for this case
    always_negative_enz_pathway =
        MetabolicPathway((:S_ext,), ((:NegativeEnzyme, (:S_ext,), (:S,)),))
    # Define an enzyme rate that incorrectly returns negative value when products are absent
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:NegativeEnzyme,(:S_ext,),(:S,)},
        metabs,
        params,
    )
        return -1.0 * params.NegativeEnzyme_Vmax * metabs.S_ext
    end
    always_negative_enz_metabs = LVector(S_ext = 1.0, S = 0.0)  # S product is absent
    always_negative_enz_params =
        LVector(NegativeEnzyme_Vmax = 1.0, NegativeEnzyme_Keq = 10.0)
    @test_throws ErrorException CellMetabolismBase.validate_enzyme_rates(
        always_negative_enz_pathway,
        always_negative_enz_metabs,
        always_negative_enz_params,
    )

    # ------- Test: Enzyme rate is negative when substrates absent, products present -------
    # Create a specialized test for this case
    always_positive_enz_pathway =
        MetabolicPathway((:S_ext,), ((:PositiveEnzyme, (:S_ext,), (:S,)),))
    # Define an enzyme rate that incorrectly returns positive value when substrates are absent
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:PositiveEnzyme,(:S_ext,),(:S,)},
        metabs,
        params,
    )
        return 1.0
    end
    always_positive_enz_metabs = LVector(S_ext = 0.0, S = 1.0)  # S_ext substrate is absent
    always_positive_enz_params =
        LVector(PositiveEnzyme_Vmax = 1.0, PositiveEnzyme_Keq = 10.0)
    # This should throw an error
    @test_throws ErrorException CellMetabolismBase.validate_enzyme_rates(
        always_positive_enz_pathway,
        always_positive_enz_metabs,
        always_positive_enz_params,
    )

    # ------- Test: Enzyme rate is zero when all substrates and products absent -------
    # Create a specialized test for this case
    never_zero_pathway = MetabolicPathway((:S_ext,), ((:NeverZeroEnz, (:S_ext,), (:S,)),))

    # Define an enzyme rate that incorrectly returns non-zero when all metabolites are absent
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:NeverZeroEnz,(:S_ext,),(:S,)},
        metabs,
        params,
    )
        # This function always returns non-zero regardless of zero concentrations
        # This violates the rule that rate should be zero when all metabolites are zero
        return params.NeverZeroEnz_Vmax
    end

    never_zero_metabs = LVector(S_ext = 0.0, S = 0.0)  # Both substrate and product absent
    never_zero_params = LVector(NeverZeroEnz_Vmax = 1.0, NeverZeroEnz_Keq = 10.0)

    # This should throw an error
    @test_throws ErrorException CellMetabolismBase.validate_enzyme_rates(
        never_zero_pathway,
        never_zero_metabs,
        never_zero_params,
    )

    # ------- Test: Enzyme rate is zero at thermodynamic equilibrium -------
    # Define a simple reversible enzyme for equilibrium testing
    non_eq_pathway = MetabolicPathway((:A_ext,), ((:NonEquilEnz, (:A_ext,), (:A,)),))
    # Define an enzyme rate that does not return zero at equilibrium
    function CellMetabolismBase.enzyme_rate(
        ::Enzyme{:NonEquilEnz,(:A_ext,),(:A,)},
        metabs,
        params,
    )
        return params.NonEquilEnz_Vmax *
               (metabs.A_ext - 1.001 * metabs.A / params.NonEquilEnz_Keq) / (
            params.NonEquilEnz_Km *
            (1 + metabs.A_ext / params.NonEquilEnz_Km + metabs.A / params.NonEquilEnz_Ki)
        )
    end

    non_eq_params = LVector(
        NonEquilEnz_Vmax = 1.0,
        NonEquilEnz_Km = 0.5,
        NonEquilEnz_Ki = 0.5,
        NonEquilEnz_Keq = 2.0,
    )

    # Create equilibrium condition: A = A_ext * Keq
    non_eq_metabs = LVector(A_ext = 1.0, A = 2.0)  # A = A_ext * Keq = 1.0 * 2.0

    # This should throw an error because function doesn't return zero at equilibrium
    @test_throws ErrorException CellMetabolismBase.validate_enzyme_rates(
        non_eq_pathway,
        non_eq_metabs,
        non_eq_params,
    )
end
