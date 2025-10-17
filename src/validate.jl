using LabelledArrays

"""
    validate(metabolic_pathway::MetabolicPathway, init_cond::LArray, params::LArray)

Perform structural checks on a metabolic pathway definition, ensuring metabolites present in the
pathway exist in the labelled initial conditions and that all enzyme rate functions behave
consistently across basic scenarios, like `rate` being positive when substrates are present
and products are absent, negative when products are present and substrates are absent, and zero
when all substrates and products are absent or at equilibrium.
"""
function validate(
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

        overlap_mid_values = Dict{Symbol,Float64}()

        for reg in regulator_list
            hasproperty(init_cond, reg) ||
                error("Regulator $(reg) for enzyme $(name(enzyme)) is missing from initial conditions; cannot validate remove_regulation.")

            test_metabs_high = @LArray eps() .+ rand(length(init_cond)) propertynames(init_cond)
            test_metabs_low = deepcopy(test_metabs_high)
            test_metabs_high[reg] = 1.0 + rand()
            test_metabs_low[reg] = 0.0

            test_params = @LArray eps() .+ rand(length(params)) propertynames(params)

            is_overlap_inhibitor =
                (reg in inhibitors(enzyme)) &&
                ((reg in substrates(enzyme)) || (reg in products(enzyme)))
            if is_overlap_inhibitor
                overlap_mid_values[reg] = _validate_inhibitor_overlap(
                    enzyme,
                    reg,
                    init_cond,
                    params,
                    test_params,
                )
                continue
            end
            if (reg in activators(enzyme)) &&
               ((reg in substrates(enzyme)) || (reg in products(enzyme)))
                overlap_mid_values[reg] = _validate_activator_overlap(
                    enzyme,
                    reg,
                    init_cond,
                    params,
                    test_params,
                )
                continue
            end

            specific_params = try
                remove_regulation(enzyme, test_params, Val(reg))
            catch err
                error("remove_regulation($(name(enzyme)), params, Val($(reg))) failed during validation: $(err)")
            end

            rate_high = rate(enzyme, test_metabs_high, specific_params)
            rate_low = rate(enzyme, test_metabs_low, specific_params)
            (isfinite(rate_high) && isfinite(rate_low)) ||
                error("remove_regulation for enzyme $(name(enzyme)) and regulator $(reg) produced non-finite rates during validation.")

            isapprox(rate_high, rate_low; atol=1e-8, rtol=1e-6) ||
                error("remove_regulation($(name(enzyme)), params, Val($(reg))) did not eliminate dependence on regulator $(reg).")
        end

        # Validate combined regulator removal
        combined_metabs_high = @LArray eps() .+ rand(length(init_cond)) propertynames(init_cond)
        combined_metabs_low = deepcopy(combined_metabs_high)
        for reg in regulator_list
            if haskey(overlap_mid_values, reg)
                combined_metabs_high[reg] = overlap_mid_values[reg]
                combined_metabs_low[reg] = overlap_mid_values[reg]
            else
                combined_metabs_high[reg] = 1.0 + rand()
                combined_metabs_low[reg] = 0.0
            end
        end
        combined_params = @LArray eps() .+ rand(length(params)) propertynames(params)

        all_removed_params = try
            remove_regulation(enzyme, combined_params)
        catch err
            error("remove_regulation($(name(enzyme)), params) failed during validation: $(err)")
        end

        rate_high_all = rate(enzyme, combined_metabs_high, all_removed_params)
        rate_low_all = rate(enzyme, combined_metabs_low, all_removed_params)
        (isfinite(rate_high_all) && isfinite(rate_low_all)) ||
            error("remove_regulation($(name(enzyme)), params) produced non-finite rates during validation.")

        isapprox(rate_high_all, rate_low_all; atol=1e-8, rtol=1e-6) ||
            error("remove_regulation($(name(enzyme)), params) did not eliminate dependence on its regulators $(Tuple(regulator_list)).")
    end

    return nothing
end

function _validate_inhibitor_overlap(
    enzyme::Enzyme,
    reg::Symbol,
    init_cond,
    params_reference,
    params_trial,
)
    params_removed = try
        remove_regulation(enzyme, deepcopy(params_trial), Val(reg))
    catch err
        error("remove_regulation($(name(enzyme)), params, Val($(reg))) failed during validation: $(err)")
    end
    _regulator_parameters_changed(params_trial, params_removed, reg) ||
        error("remove_regulation($(name(enzyme)), params, Val($(reg))) did not modify parameters associated with regulator $(reg).")
    _regulator_parameters_changed(params_trial, params_removed, reg) ||
        error("remove_regulation($(name(enzyme)), params, Val($(reg))) did not modify parameters associated with regulator $(reg).")

    mid_val = _regulator_mid_value(enzyme, params_reference, reg)
    if !(mid_val isa Real && isfinite(mid_val) && mid_val > 0)
        mid_val = 1.0
    end
    low_val = max(mid_val / 10, eps())
    high_val = max(mid_val * 10, low_val + eps())

    names = propertynames(init_cond)
    low_data = fill(mid_val, length(init_cond))
    mid_data = copy(low_data)
    high_data = copy(low_data)
    low_metabs = @LArray low_data names
    mid_metabs = @LArray mid_data names
    high_metabs = @LArray high_data names

    subs = substrates(enzyme)
    prods = products(enzyme)

    if reg in subs
        for product in prods
            low_metabs[product] = 0.0
            mid_metabs[product] = 0.0
            high_metabs[product] = 0.0
        end
        for substrate in subs
            if substrate == reg
                low_metabs[substrate] = low_val
                mid_metabs[substrate] = mid_val
                high_metabs[substrate] = high_val
            else
                low_metabs[substrate] = mid_val
                mid_metabs[substrate] = mid_val
                high_metabs[substrate] = mid_val
            end
        end
    else
        for substrate in subs
            low_metabs[substrate] = 0.0
            mid_metabs[substrate] = 0.0
            high_metabs[substrate] = 0.0
        end
        for product in prods
            if product == reg
                low_metabs[product] = low_val
                mid_metabs[product] = mid_val
                high_metabs[product] = high_val
            else
                low_metabs[product] = mid_val
                mid_metabs[product] = mid_val
                high_metabs[product] = mid_val
            end
        end
    end

    rate_low = rate(enzyme, low_metabs, params_removed)
    rate_mid = rate(enzyme, mid_metabs, params_removed)
    rate_high = rate(enzyme, high_metabs, params_removed)
    all(isfinite, (rate_low, rate_mid, rate_high)) ||
        error("remove_regulation($(name(enzyme)), params, Val($(reg))) produced non-finite rates during validation.")

    tol = 1e-8 + 1e-6 * maximum(abs, (rate_low, rate_mid, rate_high))
    if reg in subs
        ((rate_low <= rate_mid + tol) && (rate_mid <= rate_high + tol)) ||
            error("remove_regulation($(name(enzyme)), params, Val($(reg))) did not eliminate inhibitory behaviour for substrate regulator $(reg).")
    else
        ((rate_low + tol >= rate_mid) && (rate_mid + tol >= rate_high)) ||
            error("remove_regulation($(name(enzyme)), params, Val($(reg))) did not eliminate inhibitory behaviour for product regulator $(reg).")
    end

    return mid_val
end

function _validate_activator_overlap(
    enzyme::Enzyme,
    reg::Symbol,
    init_cond,
    params_reference,
    params_trial,
)
    params_removed = try
        remove_regulation(enzyme, deepcopy(params_trial), Val(reg))
    catch err
        error("remove_regulation($(name(enzyme)), params, Val($(reg))) failed during validation: $(err)")
    end

    mid_val = _regulator_mid_value(enzyme, params_reference, reg)
    if !(mid_val isa Real && isfinite(mid_val) && mid_val > 0)
        mid_val = 1.0
    end
    low_val = max(mid_val / 10, eps())
    high_val = max(mid_val * 10, low_val + eps())

    names = propertynames(init_cond)
    low_data = fill(mid_val, length(init_cond))
    high_data = copy(low_data)
    low_metabs = @LArray low_data names
    high_metabs = @LArray high_data names

    if reg in substrates(enzyme)
        for substrate in substrates(enzyme)
            if substrate == reg
                low_metabs[substrate] = low_val
                high_metabs[substrate] = high_val
            else
                low_metabs[substrate] = 0.0
                high_metabs[substrate] = 0.0
            end
        end
        for product in products(enzyme)
            low_metabs[product] = mid_val
            high_metabs[product] = mid_val
        end
    else
        for product in products(enzyme)
            if product == reg
                low_metabs[product] = low_val
                high_metabs[product] = high_val
            else
                low_metabs[product] = 0.0
                high_metabs[product] = 0.0
            end
        end
        for substrate in substrates(enzyme)
            low_metabs[substrate] = mid_val
            high_metabs[substrate] = mid_val
        end
    end

    rate_low_with = rate(enzyme, low_metabs, params_trial)
    rate_high_with = rate(enzyme, high_metabs, params_trial)
    all(isfinite, (rate_low_with, rate_high_with)) ||
        error("Rate evaluation for $(name(enzyme)) with activator $(reg) produced non-finite values during validation.")

    rate_low_removed = rate(enzyme, low_metabs, params_removed)
    rate_high_removed = rate(enzyme, high_metabs, params_removed)
    all(isfinite, (rate_low_removed, rate_high_removed)) ||
        error("remove_regulation($(name(enzyme)), params, Val($(reg))) produced non-finite rates during validation.")

    return mid_val
end

function _regulator_mid_value(enzyme::Enzyme, params, reg::Symbol)
    prefix = string(name(enzyme), "_K_")
    target = lowercase(String(reg))
    values = Float64[]
    for sym in propertynames(params)
        str = String(sym)
        startswith(str, prefix) || continue
        occursin(target, lowercase(str)) || continue
        val = params[sym]
        val isa Real || continue
        val > 0 || continue
        push!(values, float(val))
    end
    isempty(values) && return 1.0
    return exp(sum(log, values) / length(values))
end

function _regulator_parameters_changed(before, after, reg::Symbol)
    target = lowercase(String(reg))
    for sym in propertynames(before)
        hasproperty(after, sym) || continue
        occursin(target, lowercase(String(sym))) || continue
        before_val = before[sym]
        after_val = after[sym]
        if !(isequal(before_val, after_val) || (before_val isa Real && after_val isa Real && isapprox(before_val, after_val; atol=1e-12, rtol=1e-8)))
            return true
        end
    end
    return false
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
