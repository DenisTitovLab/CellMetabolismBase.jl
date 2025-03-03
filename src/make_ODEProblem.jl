using SciMLBase

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
    enz_rates = enzyme_rates(metab_path, metabs, params)
    calculate_dmetabs_from_enz_rates!(metab_path, dmetabs, enz_rates)
    return nothing
end

@inline @generated function calculate_dmetabs_from_enz_rates!(
    ::MetabolicPathway{ConstMetabs,Enzs},
    dmetabs::LArray{T,1,Vector{T},Syms},
    rates::NTuple{N},
) where {ConstMetabs,Enzs,T,Syms,N}
    temp = []
    for metab in Syms
        dmetab = :(0.0)
        if metab ∉ ConstMetabs
            for enz in Enzs
                if metab in enz[2]
                    stoich_coeff = sum(metab .== enz[2])
                    dmetab = :($dmetab - rates[$(findfirst(==(enz), Enzs))] * $stoich_coeff)
                end
                if metab in enz[3]
                    stoich_coeff = sum(metab .== enz[3])
                    dmetab = :($dmetab + rates[$(findfirst(==(enz), Enzs))] * $stoich_coeff)
                end
            end
        end
        push!(temp, :(dmetabs.$metab = $dmetab))
    end
    expr = Expr(:block, temp..., :(return nothing))
    return expr
end

function enzyme_rates(metab_path::MetabolicPathway, metabs::LArray, params::LArray)
    enzymes = _generate_Enzymes(metab_path)
    return map(enzyme -> CellMetabolismBase.rate(enzyme, metabs, params), enzymes)
end

_generate_Enzymes(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Enz -> Enzyme(Enz...), Enzs)
