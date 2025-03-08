using LabelledArrays

function validate_metabolic_pathway(
    metabolic_pathway::MetabolicPathway,
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {T1<:Real,T2<:Real,MetabNames,ParamNames}

    pathway_metab_names = Symbol[]
    #TODO: have a dedicated function to extract all metabs and also regulators
    for metabs in [substrate_names(metabolic_pathway)..., product_names(metabolic_pathway)...]
        for metab in metabs
            if !(metab in pathway_metab_names)
                push!(pathway_metab_names, metab)
            end
        end
    end
    # validate metabolic pathway metabolites are in MetabNames
    for metab in pathway_metab_names
        println("Validating metabolite: $metab")
        if !(metab in MetabNames)
            error("Metabolite $metab not found in initial conditions LArray.")
        end
    end
    validate_enzymes(metabolic_pathway, init_cond, params)

    return nothing
end

function validate_enzymes(
    metabolic_pathway::MetabolicPathway{ConstMetabs,Enzs},
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {ConstMetabs,Enzs,T1<:Real,T2<:Real,MetabNames,ParamNames}
    # Validate enzyme rate functions return correct type and check rates
    for Enz in Enzs
        rand_test_metabs = @LArray rand(length(init_cond)) propertynames(init_cond)
        enzyme_rate(Enzyme(Enz...), init_cond, params)
        substrate_names = Enz[2]
        product_names = Enz[3]
        # regulator_names = Enz[4]

        # enzyme_rate is positive if one of products is absent and substrates are present
        for product in product_names
            test_metabs = deepcopy(rand_test_metabs)
            test_metabs[product] = 0.0
            enzyme_rate(Enzyme(Enz...), test_metabs, params) > 0.0 ||
                error("Enzyme $(Enz[1]) rate should be positive when all substrates are present but products are missing.")
        end

        # enzyme_rate is negative if one of substrates is absent and products are present
        for substrate in substrate_names
            test_metabs = deepcopy(rand_test_metabs)
            test_metabs[substrate] = 0.0
            enzyme_rate(Enzyme(Enz...), test_metabs, params) < 0.0 ||
                error("Enzyme $(Enz[1]) rate should be negative when all substrates are present, indicating potential error or misconfiguration.")
        end

        # enzyme_rate is zero if all substrates and products are absent
        empty_metabs = deepcopy(rand_test_metabs)
        for metab in [substrate_names..., product_names...]
            empty_metabs[metab] = 0.0
        end
        enzyme_rate(Enzyme(Enz...), empty_metabs, params) == 0.0 || error("Enzyme $(Enz[1]) rate should be zero when all substrates and products are absent.")

        # enzyme_rate is zero when substrates and products are at equilibrium
        equilibrium_metabs = deepcopy(rand_test_metabs)
        if length(unique(product_names)) == length(product_names)
            equilibrium_metabs[product_names[1]] = params[Symbol(Enz[1], "_Keq")] * reduce(*, [equilibrium_metabs[substrate] for substrate in substrate_names], init=1.0) /
                                          reduce(*, [equilibrium_metabs[product] for product in product_names if product != product], init=1.0)
        else
            equilibrium_metabs[product_names[1]] = sqrt(params[Symbol(Enz[1], "_Keq")] * reduce(*, [equilibrium_metabs[substrate] for substrate in substrate_names], init=1.0) /
                                               reduce(*, [equilibrium_metabs[product] for product in product_names if product != product], init=1.0))
        end
        isapprox(1.0-enzyme_rate(Enzyme(Enz...), equilibrium_metabs, params), 1.0) ||
            error("Enzyme $(Enz[1]) rate should be zero when all substrates and products are at equilibrium.")
    end
    return nothing
end
