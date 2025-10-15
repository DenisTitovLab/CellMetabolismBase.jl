using SciMLBase: ODEProblem

"""
    make_ODEProblem(
    metabolic_pathway::MetabolicPathway,
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    tspan::Tuple{<:Real,<:Real},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {T1,T2,MetabNames,ParamNames}

Construct an ODEProblem for simulating a metabolic pathway.

# Arguments
- `metabolic_pathway::MetabolicPathway`: A structure that defines the metabolic pathway, including the substrates, products, and related reaction details.
- `init_cond::LArray{T,1,Vector{T},MetabNames}`: The initial conditions for the metabolites.
- `tspan::Tuple{<:Real, <:Real}`: A tuple specifying the start and end times for the simulation.
- `params::LArray{T,1,Vector{T},ParamNames}`: Parameters used in the metabolic reactions (e.g., kinetic constants).

# Returns
An `ODEProblem` instance that encapsulates the differential equations governing the metabolic pathway. This problem can be solved using ODE solvers from SciMLBase.

# Details
The returned ODEProblem uses `metabolicpathway_odes!` to compute the time evolution of metabolite concentrations
by adjusting their rates based on enzyme kinetics defined in `enzyme_rates`. Ensure that `metabolic_pathway`,
`init_cond`, and `params` have matching entries for substrates, products, and parameters for correct simulation.

# Example
```julia
using CellMetabolismBase
using DifferentialEquations

metab_path = MetabolicPathway(...)
init_cond = LArray(...)
tspan = (0.0, 1e8)
params = LArray(...)
prob = make_ODEProblem(metab_path, init_cond, tspan, params)
sol = solve(prob, RadauIIA9(), abstol=1e-15, reltol=1e-8)
```
"""
function make_ODEProblem(
    metabolic_pathway::MetabolicPathway,
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    tspan::Tuple{<:Real,<:Real},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {T1<:Real,T2<:Real,MetabNames,ParamNames}
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
    dmetabs::LArray{T1,1,Vector{T1},MetabNames},
    metabs::LArray{T2,1,Vector{T2},MetabNames},
    params::LArray{T3,1,Vector{T3},ParamNames},
    t,
) where {T1<:Real,T2<:Real,T3<:Real,MetabNames,ParamNames}
    enz_rates = enzyme_rates(metab_path, metabs, params)
    calculate_dmetabs_from_enz_rates!(metab_path, dmetabs, enz_rates)
    return nothing
end

@inline @generated function calculate_dmetabs_from_enz_rates!(
    ::MetabolicPathway{ConstMetabs,Enzs},
    dmetabs::LArray{T,1,Vector{T},MetabNames},
    rates::NTuple{N},
) where {ConstMetabs,Enzs,T<:Real,MetabNames,N}
    temp = []
    for metab in MetabNames
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

function enzyme_rates(
    metab_path::MetabolicPathway,
    metabs::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {T1<:Real,T2<:Real,MetabNames,ParamNames}
    enzymes = _generate_Enzymes(metab_path)
    return map(enzyme -> CellMetabolismBase.enzyme_rate(enzyme, metabs, params), enzymes)
end

"""
    disequilibrium_ratios(metab_path::MetabolicPathway,
        metabs::LArray{T1,1,Vector{T1},MetabNames},
        params::LArray{T2,1,Vector{T2},ParamNames})

Return the disequilibrium (reaction quotient over equilibrium constant) for each enzyme
in `metab_path`, ordered according to `enzyme_names(metab_path)`.
"""
@generated function disequilibrium_ratios(
    ::MetabolicPathway{ConstMetabs,Enzs},
    metabs::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {ConstMetabs,Enzs,T1<:Real,T2<:Real,MetabNames,ParamNames}
    ratio_exprs = Expr[]
    for (idx, Enz) in enumerate(Enzs)
        enz_name = Enz[1]
        substrates = Enz[2]
        products = Enz[3]
        keq_sym = Symbol(enz_name, "_Keq")
        error_msg = "Parameter $(keq_sym) (equilibrium constant) not found for enzyme $(enz_name). Please ensure the parameter is defined in the params object."

        numerator_expr = :(one(eltype(metabs)))
        for product in products
            numerator_expr = :( $numerator_expr * metabs.$product )
        end

        denominator_expr = :(one(eltype(metabs)))
        for substrate in substrates
            denominator_expr = :( $denominator_expr * metabs.$substrate )
        end

        reaction_quotient_sym = gensym(:reaction_quotient)
        block_expr = Expr(
            :block,
            :(hasproperty(params, $(QuoteNode(keq_sym))) ||
                throw(ArgumentError($(QuoteNode(error_msg))))),
            Expr(:(=), reaction_quotient_sym, :( $numerator_expr / $denominator_expr )),
            :( $reaction_quotient_sym / getproperty(params, $(QuoteNode(keq_sym))) ),
        )
        push!(ratio_exprs, block_expr)
    end
    return Expr(:tuple, ratio_exprs...)
end

@generated _generate_Enzymes(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Enz -> Enzyme(Enz...), Enzs)
