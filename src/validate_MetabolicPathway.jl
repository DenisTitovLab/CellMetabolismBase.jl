using LabelledArrays

"""
    validate_MetabolicPathway(metabolic_pathway::MetabolicPathway, init_cond::LArray, params::LArray)

Perform structural checks on a metabolic pathway definition, ensuring metabolites present in the
pathway exist in the labelled initial conditions and that all enzyme rate functions behave
consistently across basic scenarios, like `rate` being positive when substrates are present
and products are absent, negative when products are present and substrates are absent, and zero
when all substrates and products are absent or at equilibrium.
"""
function validate_MetabolicPathway(
    metabolic_pathway::MetabolicPathway,
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {T1<:Real,T2<:Real,MetabNames,ParamNames}

    # validate metabolic pathway metabolites are in MetabNames
    for metab in metabolites(metabolic_pathway)
        if !(metab in MetabNames)
            error("Metabolite $metab not found in initial conditions LArray.")
        end
    end

    #=
    TODO: validate parameter names
        - enzymes only use params that start with enzyme name like :Enz_K_a_S1
        - maybe enforce that params also have Metabolite names in params that correspond to enzyme substrates/products/regulators
    =#

    #validate rate equations
    validate_enzyme_rates(metabolic_pathway, init_cond, params)

    #validate regulation removal behaviour when regulators are present
    validate_regulation_removal(metabolic_pathway, init_cond, params)

    return nothing
end

function validate_regulation_removal(
    metabolic_pathway::MetabolicPathway,
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {T1<:Real,T2<:Real,MetabNames,ParamNames}
    enzymes = CellMetabolismBase._generate_Enzymes(metabolic_pathway)
    isempty(enzymes) && return nothing

    for enzyme in enzymes
        regulator_list = Symbol[]
        append!(regulator_list, activators(enzyme))
        append!(regulator_list, inhibitors(enzyme))
        isempty(regulator_list) && continue

        for reg in regulator_list
            hasproperty(init_cond, reg) ||
                error("Regulator $(reg) for enzyme $(name(enzyme)) is missing from initial conditions; cannot validate remove_regulation.")

            test_metabs_high = @LArray eps() .+ rand(length(init_cond)) propertynames(init_cond)
            test_metabs_low = deepcopy(test_metabs_high)
            test_metabs_high[reg] = 1.0 + rand()
            test_metabs_low[reg] = 0.0

            test_params = @LArray eps() .+ rand(length(params)) propertynames(params)

            specific_params = try
                remove_regulation(test_params, enzyme, Val(reg))
            catch err
                error("remove_regulation(params, $(name(enzyme)), Val($(reg))) failed during validation: $(err)")
            end

            rate_high = rate(enzyme, test_metabs_high, specific_params)
            rate_low = rate(enzyme, test_metabs_low, specific_params)
            (isfinite(rate_high) && isfinite(rate_low)) ||
                error("remove_regulation for enzyme $(name(enzyme)) and regulator $(reg) produced non-finite rates during validation.")

            isapprox(rate_high, rate_low; atol=1e-8, rtol=1e-6) ||
                error("remove_regulation(params, $(name(enzyme)), Val($(reg))) did not eliminate dependence on regulator $(reg).")
        end

        # Validate combined regulator removal
        combined_metabs_high = @LArray eps() .+ rand(length(init_cond)) propertynames(init_cond)
        combined_metabs_low = deepcopy(combined_metabs_high)
        for reg in regulator_list
            combined_metabs_high[reg] = 1.0 + rand()
            combined_metabs_low[reg] = 0.0
        end
        combined_params = @LArray eps() .+ rand(length(params)) propertynames(params)

        all_removed_params = try
            remove_regulation(combined_params, enzyme)
        catch err
            error("remove_regulation(params, $(name(enzyme))) failed during validation: $(err)")
        end

        rate_high_all = rate(enzyme, combined_metabs_high, all_removed_params)
        rate_low_all = rate(enzyme, combined_metabs_low, all_removed_params)
        (isfinite(rate_high_all) && isfinite(rate_low_all)) ||
            error("remove_regulation(params, $(name(enzyme))) produced non-finite rates during validation.")

        isapprox(rate_high_all, rate_low_all; atol=1e-8, rtol=1e-6) ||
            error("remove_regulation(params, $(name(enzyme))) did not eliminate dependence on its regulators $(Tuple(regulator_list)).")
    end

    return nothing
end

function validate_enzyme_rates(
    ::MetabolicPathway{ConstMetabs,Enzs},
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {ConstMetabs,Enzs,T1<:Real,T2<:Real,MetabNames,ParamNames}
    # Validate enzyme rate functions return correct type and check rates
    for Enz in Enzs
        rand_test_metabs = @LArray eps() .+ rand(length(init_cond)) propertynames(init_cond)
        rand_params = @LArray eps() .+ rand(length(params)) propertynames(params)
        test_rate = rate(Enzyme(Enz...), init_cond, rand_params)
        typeof(test_rate) <: Real ||
            error("Enzyme $(Enz[1]) rate function should return a Real number.")
        substrate_names = Enz[2]
        product_names = Enz[3]
        # regulator_names = Enz[4]

        # rate is positive if one of products is absent and substrates are present
        for product in product_names
            test_metabs = deepcopy(rand_test_metabs)
            test_metabs[product] = 0.0
            rate(Enzyme(Enz...), test_metabs, rand_params) > 0.0 ||
                error("Enzyme $(Enz[1]) rate should be positive when substrates are present and one product is missing.")
        end

        # rate is negative if one of substrates is absent and products are present
        for substrate in substrate_names
            test_metabs = deepcopy(rand_test_metabs)
            test_metabs[substrate] = 0.0
            rate(Enzyme(Enz...), test_metabs, rand_params) < 0.0 ||
                error("Enzyme $(Enz[1]) rate should be negative when products are present and one substrate is missing.")
        end

        # rate is zero if all substrates and products are absent
        empty_metabs = deepcopy(rand_test_metabs)
        for metab in [substrate_names..., product_names...]
            empty_metabs[metab] = 0.0
        end
        rate(Enzyme(Enz...), empty_metabs, rand_params) == 0.0 || error("Enzyme $(Enz[1]) rate should be zero when all substrates and products are absent.")

        # rate is zero when substrates and products are at equilibrium
        equilibrium_metabs = deepcopy(rand_test_metabs)
        if length(unique(product_names)) == length(product_names)
            equilibrium_metabs[product_names[1]] = rand_params[Symbol(Enz[1], "_Keq")] * reduce(*, [equilibrium_metabs[substrate] for substrate in substrate_names], init=1.0) /
                                          reduce(*, [equilibrium_metabs[product] for product in product_names[2:end]], init=1.0)
        else
            equilibrium_metabs[product_names[1]] = sqrt(rand_params[Symbol(Enz[1], "_Keq")] * reduce(*, [equilibrium_metabs[substrate] for substrate in substrate_names], init=1.0) /
                                               reduce(*, [equilibrium_metabs[product] for product in product_names if product != product], init=1.0))
        end
        isapprox(1.0 - rate(Enzyme(Enz...), equilibrium_metabs, rand_params), 1.0) ||
            error("Enzyme $(Enz[1]) rate = $(rate(Enzyme(Enz...), equilibrium_metabs, rand_params)) when substrates and products are at equilibrium but should be zero.")
    end
    return nothing
end
