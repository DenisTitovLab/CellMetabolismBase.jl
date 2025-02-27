#TODO: add dependence on SciMLBase or wherever ODEProblem is defined
using SciMLBase
# include("types_and_related_functions.jl")



"""
    make_ODEProblem(metabolic_pathway, init_cond, tspan, params)

Construct an ODEProblem for simulating a metabolic pathway.

# Arguments
- `metabolic_pathway::MetabolicPathway`: A structure that defines the metabolic pathway, including the substrates, products, and related reaction details.
- `init_cond::LArray`: The initial conditions for the metabolites.
- `tspan::Tuple{<:Number, <:Number}`: A tuple specifying the start and end times for the simulation.
- `params::LArray`: Parameters used in the metabolic reactions (e.g., kinetic constants).

# Returns
An `ODEProblem` instance that encapsulates the differential equations governing the metabolic pathway. This problem can be solved using ODE solvers from SciMLBase.

# Details
The returned ODEProblem uses `metabolicpathway_odes!` to compute the time evolution of metabolite concentrations
by adjusting their rates based on enzyme kinetics defined in `enzyme_rates`. Ensure that `metabolic_pathway`,
`init_cond`, and `params` have matching entries for substrates, products, and parameters for correct simulation.

"""
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
