#TODO: add dependence on SciMLBase or wherever ODEProblem is defined
using SciMLBase
include("types_and_related_functions.jl")

function make_ODEProblem(metabolic_pathway, init_cond, tspan, params)
    #test that the pathway was assembled correctly
    #assert that metabs names overlap with metabolic_pathway substrates and products
    #assert that params overlap with metabolic_pathway params or maybe use the latter?

    return ODEProblem(
        (dm, m, p, t) -> metabolicpathway_odes!(metabolic_pathway, dm, m, p, t),
        init_cond,
        tspan,
        params,
    )
end

function metabolicpathway_odes!(
    metab_path::MetabolicPathway,
    dmetabs::LArray,
    metabs::LArray,
    params::LArray,
    t,
)
    propertynames(metabs) == propertynames(dmetabs) ||
        error("metabs and dmetabs must have the same propertynames")
    dmetabs .= zero(eltype(dmetabs))
    substr_indxs = _substrate_indxs(metabs, metab_path)
    prod_indxs = _product_indxs(metabs, metab_path)
    rates = enzyme_rates(metab_path, metabs, params)
    @inbounds for (e, rate) in enumerate(rates)
        for i in substr_indxs[e]
            dmetabs[i] -= rate
        end
        for k in prod_indxs[e]
            dmetabs[k] += rate
        end
    end
    return nothing
end

function enzyme_rates(metab_path::MetabolicPathway, metabs::LArray, params::LArray)
    enzymes = _generate_enzymes(metab_path)
    return map(enzyme -> rate(enzyme, metabs, params), enzymes)
end
