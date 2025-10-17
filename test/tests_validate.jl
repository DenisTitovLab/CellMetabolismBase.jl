@testitem "MetabolicPathway Validation Tests" setup = [TestMetabolicPathway] begin
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
    @test_nowarn CellMetabolismBase.validate(test_pathway, valid_metabs, valid_params)
    # Test the enzyme validation function works
    @test_nowarn CellMetabolismBase.validate_enzyme_rates(
        test_pathway,
        valid_metabs,
        valid_params,
    )

    # ------- Test: Check metabolites found in initial conditions -------
    # Test 1: Missing metabolite - should fail
    missing_metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0) # Missing C and D
    @test_throws ErrorException CellMetabolismBase.validate(
        test_pathway,
        missing_metabs,
        valid_params,
    )

    # Test 2: Extra metabolite - should pass
    extra_metabs = LVector(A_media = 2.0, A = 1.0, B = 1.0, C = 1.0, D = 1.0, E = 0.5) # Extra E
    @test_nowarn CellMetabolismBase.validate(test_pathway, extra_metabs, valid_params)

    # ------- Test: rate() is defined and returns Real number -------
    # Test with a simple pathway
    simple_pathway = MetabolicPathway((:A_ext,), ((:SimpleEnz, (:A_ext,), (:B,)),))
    function CellMetabolismBase.rate(::Enzyme{:SimpleEnz,(:A_ext,),(:B,)}, metabs, params)
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
    function CellMetabolismBase.rate(
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
    function CellMetabolismBase.rate(
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
    function CellMetabolismBase.rate(
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

    # Should throw an error when rate function is not defined
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
    function CellMetabolismBase.rate(
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
    function CellMetabolismBase.rate(
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

    # ------- Test: Regulation removal validation -------
    regulated_pathway = MetabolicPathway((), ((:RegEnz, (:S,), (:P,), (:Act,), (:Inh,)),))

    function CellMetabolismBase.rate(
        ::Enzyme{:RegEnz,(:S,),(:P,),(:Act,),(:Inh,)},
        metabs,
        params,
    )
        driving = params.RegEnz_Vmax * (metabs.S - metabs.P / params.RegEnz_Keq)
        activation = 1 + params.RegEnz_K_act * metabs.Act
        inhibition = 1 / (1 + params.RegEnz_K_inh * metabs.Inh)
        return driving * activation * inhibition
    end

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:RegEnz,(:S,),(:P,),(:Act,),(:Inh,)},
        params,
        ::Val{:Act},
    )
        params = deepcopy(params)
        setproperty!(params, :RegEnz_K_act, 0.0)
        return params
    end

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:RegEnz,(:S,),(:P,),(:Act,),(:Inh,)},
        params,
        ::Val{:Inh},
    )
        params = deepcopy(params)
        setproperty!(params, :RegEnz_K_inh, 0.0)
        return params
    end

    regulated_metabs = LVector(S = 1.0, P = 0.2, Act = 0.4, Inh = 0.3)
    regulated_params =
        LVector(RegEnz_Vmax = 2.0, RegEnz_Keq = 5.0, RegEnz_K_act = 1.2, RegEnz_K_inh = 0.8)

    @test_throws ErrorException CellMetabolismBase.validate_regulation_removal(
        regulated_pathway,
        regulated_metabs,
        regulated_params,
    )

    function CellMetabolismBase.remove_regulation(
        enzyme::Enzyme{:RegEnz,(:S,),(:P,),(:Act,),(:Inh,)},
        params,
    )
        params = deepcopy(params)
        params = CellMetabolismBase.remove_regulation(enzyme, params, Val(:Act))
        params = CellMetabolismBase.remove_regulation(enzyme, params, Val(:Inh))
        return params
    end

    @test_nowarn CellMetabolismBase.validate(
        regulated_pathway,
        regulated_metabs,
        regulated_params,
    )

    faulty_pathway = MetabolicPathway((), ((:FaultyEnz, (:S,), (:P,), (:Act,), ()),))

    function CellMetabolismBase.rate(
        ::Enzyme{:FaultyEnz,(:S,),(:P,),(:Act,),()},
        metabs,
        params,
    )
        return params.FaultyEnz_Vmax *
               (1 + params.FaultyEnz_K_act * metabs.Act) *
               (metabs.S - metabs.P / params.FaultyEnz_Keq)
    end

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:FaultyEnz,(:S,),(:P,),(:Act,),()},
        params,
        ::Val{:Act},
    )
        return deepcopy(params) # regulator is not actually removed
    end

    function CellMetabolismBase.remove_regulation(
        enzyme::Enzyme{:FaultyEnz,(:S,),(:P,),(:Act,),()},
        params,
    )
        params = deepcopy(params)
        params = CellMetabolismBase.remove_regulation(enzyme, params, Val(:Act))
        return params
    end

    faulty_metabs = LVector(S = 1.0, P = 0.2, Act = 0.5)
    faulty_params =
        LVector(FaultyEnz_Vmax = 1.5, FaultyEnz_Keq = 3.0, FaultyEnz_K_act = 0.7)

    @test_throws ErrorException CellMetabolismBase.validate(
        faulty_pathway,
        faulty_metabs,
        faulty_params,
    )

    missing_specific_pathway =
        MetabolicPathway((), ((:MissingSpecific, (:S,), (:P,), (:Act,), ()),))

    function CellMetabolismBase.rate(
        ::Enzyme{:MissingSpecific,(:S,),(:P,),(:Act,),()},
        metabs,
        params,
    )
        return params.MissingSpecific_Vmax *
               (metabs.S - metabs.P / params.MissingSpecific_Keq) *
               (1 + params.MissingSpecific_K_act * metabs.Act)
    end

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:MissingSpecific,(:S,),(:P,),(:Act,),()},
        params,
    )
        params = deepcopy(params)
        setproperty!(params, :MissingSpecific_K_act, 0.0)
        return params
    end

    missing_specific_metabs = LVector(S = 1.0, P = 0.1, Act = 0.6)
    missing_specific_params = LVector(
        MissingSpecific_Vmax = 2.0,
        MissingSpecific_Keq = 4.0,
        MissingSpecific_K_act = 1.1,
    )

    @test_throws ErrorException CellMetabolismBase.validate(
        missing_specific_pathway,
        missing_specific_metabs,
        missing_specific_params,
    )

    # ------- Test: inhibitor is also a substrate -------
    substrate_inhibitor_pathway =
        MetabolicPathway((), ((:SubInhibEnz, (:S,), (:P,), (), (:S,)),))

    function CellMetabolismBase.rate(
        ::Enzyme{:SubInhibEnz,(:S,),(:P,),(),(:S,)},
        metabs,
        params,
    )
        driving = params.SubInhibEnz_Vmax * (metabs.S - metabs.P / params.SubInhibEnz_Keq)
        binding = 1 / (1 + metabs.S / params.SubInhibEnz_K_S)
        inhibition = 1 / (1 + metabs.S / params.SubInhibEnz_K_S_inh)
        return driving * binding * inhibition
    end

    substrate_inhibitor_metabs = LVector(S = 0.5, P = 0.0)
    substrate_inhibitor_params = LVector(
        SubInhibEnz_Vmax = 1.0,
        SubInhibEnz_Keq = 5.0,
        SubInhibEnz_K_S = 0.3,
        SubInhibEnz_K_S_inh = 0.2,
    )

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:SubInhibEnz,(:S,),(:P,),(),(:S,)},
        params,
        ::Val{:S},
    )
        return deepcopy(params)
    end

    @test_throws ErrorException CellMetabolismBase.validate(
        substrate_inhibitor_pathway,
        substrate_inhibitor_metabs,
        substrate_inhibitor_params,
    )

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:SubInhibEnz,(:S,),(:P,),(),(:S,)},
        params,
        ::Val{:S},
    )
        params = deepcopy(params)
        setproperty!(params, :SubInhibEnz_K_S_inh, Inf)
        return params
    end

    function CellMetabolismBase.remove_regulation(
        enzyme::Enzyme{:SubInhibEnz,(:S,),(:P,),(),(:S,)},
        params,
    )
        return CellMetabolismBase.remove_regulation(enzyme, deepcopy(params), Val(:S))
    end

    @test_nowarn CellMetabolismBase.validate(
        substrate_inhibitor_pathway,
        substrate_inhibitor_metabs,
        substrate_inhibitor_params,
    )

    # ------- Test: inhibitor is also a product -------
    product_inhibitor_pathway =
        MetabolicPathway((), ((:ProdInhibEnz, (:S,), (:P,), (), (:P,)),))

    function CellMetabolismBase.rate(
        ::Enzyme{:ProdInhibEnz,(:S,),(:P,),(),(:P,)},
        metabs,
        params,
    )
        driving = params.ProdInhibEnz_Vmax * (metabs.S - metabs.P / params.ProdInhibEnz_Keq)
        binding =
            1 /
            (1 + metabs.S / params.ProdInhibEnz_K_S + metabs.P / params.ProdInhibEnz_K_P)
        inhibition = 1 / (1 + metabs.P / params.ProdInhibEnz_K_P_inh)
        return driving * binding * inhibition
    end

    product_inhibitor_metabs = LVector(S = 1.0, P = 0.5)
    product_inhibitor_params = LVector(
        ProdInhibEnz_Vmax = 1.0,
        ProdInhibEnz_Keq = 2.0,
        ProdInhibEnz_K_S = 0.4,
        ProdInhibEnz_K_P = 0.6,
        ProdInhibEnz_K_P_inh = 0.3,
    )

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:ProdInhibEnz,(:S,),(:P,),(),(:P,)},
        params,
        ::Val{:P},
    )
        return deepcopy(params)
    end

    @test_throws ErrorException CellMetabolismBase.validate(
        product_inhibitor_pathway,
        product_inhibitor_metabs,
        product_inhibitor_params,
    )

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:ProdInhibEnz,(:S,),(:P,),(),(:P,)},
        params,
        ::Val{:P},
    )
        params = deepcopy(params)
        setproperty!(params, :ProdInhibEnz_K_P_inh, Inf)
        return params
    end

    function CellMetabolismBase.remove_regulation(
        enzyme::Enzyme{:ProdInhibEnz,(:S,),(:P,),(),(:P,)},
        params,
    )
        return CellMetabolismBase.remove_regulation(enzyme, deepcopy(params), Val(:P))
    end

    @test_nowarn CellMetabolismBase.validate(
        product_inhibitor_pathway,
        product_inhibitor_metabs,
        product_inhibitor_params,
    )

    # ------- Test: activator is also a substrate -------
    substrate_activator_pathway =
        MetabolicPathway((), ((:SubActEnz, (:S,), (:P,), (:S,), ()),))

    function CellMetabolismBase.rate(
        ::Enzyme{:SubActEnz,(:S,),(:P,),(:S,),()},
        metabs,
        params,
    )
        driving = params.SubActEnz_Vmax * (metabs.S - metabs.P / params.SubActEnz_Keq)
        activation = 1 + params.SubActEnz_K_act * metabs.S
        return driving * activation
    end

    substrate_activator_metabs = LVector(S = 0.4, P = 1.2)
    substrate_activator_params =
        LVector(SubActEnz_Vmax = 1.0, SubActEnz_Keq = 0.2, SubActEnz_K_act = 2.0)

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:SubActEnz,(:S,),(:P,),(:S,),()},
        params,
        ::Val{:S},
    )
        return deepcopy(params)
    end

    @test_throws ErrorException CellMetabolismBase.validate(
        substrate_activator_pathway,
        substrate_activator_metabs,
        substrate_activator_params,
    )

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:SubActEnz,(:S,),(:P,),(:S,),()},
        params,
        ::Val{:S},
    )
        params = deepcopy(params)
        setproperty!(params, :SubActEnz_K_act, 0.0)
        return params
    end

    function CellMetabolismBase.remove_regulation(
        enzyme::Enzyme{:SubActEnz,(:S,),(:P,),(:S,),()},
        params,
    )
        return CellMetabolismBase.remove_regulation(enzyme, deepcopy(params), Val(:S))
    end

    @test_nowarn CellMetabolismBase.validate(
        substrate_activator_pathway,
        substrate_activator_metabs,
        substrate_activator_params,
    )

    # ------- Test: activator is also a product -------
    product_activator_pathway =
        MetabolicPathway((), ((:ProdActEnz, (:S,), (:P,), (:P,), ()),))

    function CellMetabolismBase.rate(
        ::Enzyme{:ProdActEnz,(:S,),(:P,),(:P,),()},
        metabs,
        params,
    )
        driving = params.ProdActEnz_Vmax * (metabs.S - metabs.P / params.ProdActEnz_Keq)
        activation = 1 + params.ProdActEnz_K_act * metabs.P
        return driving * activation
    end

    product_activator_metabs = LVector(S = 1.0, P = 0.2)
    product_activator_params =
        LVector(ProdActEnz_Vmax = 1.2, ProdActEnz_Keq = 5.0, ProdActEnz_K_act = 3.0)

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:ProdActEnz,(:S,),(:P,),(:P,),()},
        params,
        ::Val{:P},
    )
        return deepcopy(params)
    end

    @test_throws ErrorException CellMetabolismBase.validate(
        product_activator_pathway,
        product_activator_metabs,
        product_activator_params,
    )

    function CellMetabolismBase.remove_regulation(
        ::Enzyme{:ProdActEnz,(:S,),(:P,),(:P,),()},
        params,
        ::Val{:P},
    )
        params = deepcopy(params)
        setproperty!(params, :ProdActEnz_K_act, 0.0)
        return params
    end

    function CellMetabolismBase.remove_regulation(
        enzyme::Enzyme{:ProdActEnz,(:S,),(:P,),(:P,),()},
        params,
    )
        return CellMetabolismBase.remove_regulation(enzyme, deepcopy(params), Val(:P))
    end

    @test_nowarn CellMetabolismBase.validate(
        product_activator_pathway,
        product_activator_metabs,
        product_activator_params,
    )

    # ------- Test: Enzyme rate is zero when all substrates and products absent -------
    # Create a specialized test for this case
    never_zero_pathway = MetabolicPathway((:S_ext,), ((:NeverZeroEnz, (:S_ext,), (:S,)),))

    # Define an enzyme rate that incorrectly returns non-zero when all metabolites are absent
    function CellMetabolismBase.rate(
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
    function CellMetabolismBase.rate(::Enzyme{:NonEquilEnz,(:A_ext,),(:A,)}, metabs, params)
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
