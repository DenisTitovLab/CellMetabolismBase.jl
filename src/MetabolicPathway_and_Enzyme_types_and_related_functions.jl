using LabelledArrays

"""
    MetabolicPathway{ConstantMetabolites,Enzymes}

Internal parametric type backing [`MetabolicPathway`](@ref). Prefer the constructor helper when
creating pathways.
"""
struct MetabolicPathway{ConstantMetabolites,Enzymes} end

#TODO: transition to using DSL and move all the synthax checking to the DSL
"""
    MetabolicPathway(constant_metabs::Tuple{Vararg{Symbol}},
                     enzymes::Tuple{Vararg{<:Tuple}})

Construct an immutable metabolic pathway definition suitable for downstream helpers such as
[`make_ODEProblem`](@ref) and [`make_EnsembleProblem`](@ref).

# Arguments
- `constant_metabs`: tuple of metabolite symbols treated as constant during simulations.
- `enzymes`: tuple of enzyme specification tuples accepted by [`Enzyme`](@ref). Each element must be
  either `(name, substrates, products)` or `(name, substrates, products, activators, inhibitors)`.

# Returns
- `MetabolicPathway{ConstantMetabolites,Enzymes}()`: a canonical pathway instance storing metabolite
  and enzyme names at the type level.

# Example
```julia
const_metabs = (:Glucose_media,)
enzyme_specs = (
    (:GLUT, (:Glucose_media,), (:Glucose,)),
    (:HK1, (:Glucose, :ATP), (:G6P, :ADP), (:Phosphate,), (:G6P,)),
)
pathway = MetabolicPathway(const_metabs, enzyme_specs)
```
"""
function MetabolicPathway(ConstantMetabolites, Enzymes)
    ConstantMetabolites isa Tuple{Vararg{Symbol}} || error("ConstantMetabolites must be a tuple of symbols like (:Glucose, :Lactate, :ATP,)")
    for e in Enzymes
        (e isa Tuple{
            Symbol,
            Tuple{Symbol,Vararg{Symbol}},
            Tuple{Symbol,Vararg{Symbol}},
            Tuple{Vararg{Symbol}},
            Tuple{Vararg{Symbol}}
        } ||
         e isa Tuple{
            Symbol,
            Tuple{Symbol,Vararg{Symbol}},
            Tuple{Symbol,Vararg{Symbol}}
        }) ||
            error("Enzymes must be a tuple of tuples like ((:Name, (:S1,), (:P1,:P2),(:A1,),(:I,)), (:Name2, (:Substrate3,), (:Product4,)),...)")
    end
    return MetabolicPathway{ConstantMetabolites,Enzymes}()
end

"""
    Enzyme{Name,Substrates,Products,Activators,Inhibitors}

Internal parametric type backing [`Enzyme`](@ref). Use the exported constructor to create instances.
"""
struct Enzyme{Name,Substrates,Products,Activators,Inhibitors} end

#TODO: transition to using DSL and move all the synthax checking to the DSL
"""
    Enzyme(name::Symbol,
           substrates::Tuple{Symbol,Vararg{Symbol}},
           products::Tuple{Symbol,Vararg{Symbol}})

    Enzyme(name::Symbol,
           substrates::Tuple{Symbol,Vararg{Symbol}},
           products::Tuple{Symbol,Vararg{Symbol}},
           activators::Tuple{Vararg{Symbol}},
           inhibitors::Tuple{Vararg{Symbol}})

Validate enzyme metadata before constructing an `Enzyme`. Ensures the presence of at least one
substrate and product while allowing optional activator and inhibitor tuples.

# Arguments
- `name`: unique symbol identifying the enzyme.
- `substrates`: tuple of substrate metabolite symbols; must contain at least one entry.
- `products`: tuple of product metabolite symbols; must contain at least one entry.
- `activators`: (optional) tuple of activator symbols, defaults to `()`.
- `inhibitors`: (optional) tuple of inhibitor symbols, defaults to `()`.

# Returns
- `Enzyme{Name,Substrates,Products,Activators,Inhibitors}()`: immutable enzyme descriptor with names
  stored at the type level.

# Example
```julia
Enzyme(:GLUT, (:Glucose_media,), (:Glucose,))
Enzyme(:HK1, (:Glucose, :ATP), (:G6P, :ADP), (:Phosphate,), (:G6P,))
```
"""
function Enzyme(Name, Substrates, Products, Activators, Inhibitors)
    Name isa Symbol || error("Name must be a symbol like :Enz1")
    Substrates isa Tuple{Symbol,Vararg{Symbol}} || error("Substrates must be a tuple of symbols like (:S1,) and enzymes must have at least one substrate")
    Products isa Tuple{Symbol,Vararg{Symbol}} || error("Products must be a tuple of symbols like (:P1,) and enzymes must have at least one product")
    Activators isa Tuple{Vararg{Symbol}} || error("Activators must be a tuple of symbols like (:A1,) or empty tuple ()")
    Inhibitors isa Tuple{Vararg{Symbol}} || error("Inhibitors must be a tuple of symbols like (:I1,) or empty tuple ()")
    return Enzyme{Name,Substrates,Products,Activators,Inhibitors}()
end
Enzyme(Name, Substrates, Products) = Enzyme(Name, Substrates, Products, (), ())
"""
    rate(enzyme::Enzyme, metabs, params)

Fallback method for enzyme rate calculations. Packages depending on `CellMetabolismBase` should extend
this function for relevant enzyme variants. The default implementation throws an error to highlight
missing rate equations.
"""
rate(enzyme::Enzyme, metabs, params) = error("rate function not defined for enzyme: $enzyme")

"""
    remove_regulation(enzyme::Enzyme, params)

Fallback method for removing all regulation from an enzyme. Packages depending on `CellMetabolismBase`
should extend this function for relevant enzyme variants. The default implementation throws an error
to highlight missing implementations.
"""
remove_regulation(enzyme::Enzyme, params) =
    error("remove_regulation not defined for enzyme: $enzyme")

"""
    remove_regulation(enzyme::Enzyme, params, ::Val{regulator})

Fallback method for removing specific regulator from an enzyme. Packages depending on `CellMetabolismBase`
should extend this function for relevant enzyme-regulator combinations. The default implementation throws
an error to highlight missing implementations.
"""
remove_regulation(enzyme::Enzyme, params, ::Val{Reg}) where {Reg} =
    error("remove_regulation not defined for enzyme: $enzyme, regulator: $Reg")

"""
    name(enzyme::Enzyme)

Return `Symbol` of enzyme name.
"""
@inline name(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Name

"""
    substrates(enzyme::Enzyme)

Return the tuple of `Symbols` of substrate names of the enzyme reaction.
"""
@inline substrates(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Substrates

"""
    products(enzyme::Enzyme)

Return the tuple of `Symbols` of product names of the enzyme reaction.
"""
@inline products(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Products

"""
    activators(enzyme::Enzyme)

Return the tuple of `Symbols` of activator names of the enzyme.
"""
@inline activators(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Activators

"""
    inhibitors(enzyme::Enzyme)

Return the tuple of `Symbols` of inhibitor names of the enzyme.
"""
@inline inhibitors(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Inhibitors

"""
    disequilibrium_ratio(enzyme::Enzyme, metabs::LArray, params::LArray)

Compute the ratio of the product of product concentrations to the product of substrate concentrations, scaled by the equilibrium constant (i.e., disequilibrium ratio) for the enzyme. Throws an `ArgumentError` if the corresponding `_Keq` parameter is missing in `params`.
"""
function disequilibrium_ratio(
    enzyme::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
    metabs::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {Name,Substrates,Products,Activators,Inhibitors,T1<:Real,T2<:Real,MetabNames,ParamNames}
    keq_sym = Symbol(Name, "_Keq")
    error_msg = "Parameter $(keq_sym) (equilibrium constant) not found for enzyme $(Name). Please ensure the parameter is defined in the params object."
    hasproperty(params, keq_sym) ||
        throw(ArgumentError(error_msg))

    numerator = one(eltype(metabs))
    for product in Products
        numerator *= getproperty(metabs, product)
    end

    denominator = one(eltype(metabs))
    for substrate in Substrates
        denominator *= getproperty(metabs, substrate)
    end

    return numerator / denominator / getproperty(params, keq_sym)
end

@generated function _generate_Enzymes(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs}
    return map(Enz -> Enzyme(Enz...), Enzs)
end

"""
    constant_metabolites(pathway::MetabolicPathway)

Return the tuple of `Symbol` of metabolite names treated as constants in the pathway definition.
"""
constant_metabolites(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = ConstMetabs

"""
    enzymes(pathway::MetabolicPathway)

Return an `NTuple` with `Symbols` of enzyme names in pathway order.
"""
function enzymes(
    pathway::MetabolicPathway,
)
    enzyme_instances = _generate_Enzymes(pathway)
    return map(name, enzyme_instances)
end

"""
    substrates(pathway::MetabolicPathway)

Return an `NTuple` where each entry holds the tuple of `Symbols` of substrate names for each enzyme, matching
the ordering produced by `enzymes(pathway)`.
"""
function substrates(
    pathway::MetabolicPathway,
)
    enzyme_instances = _generate_Enzymes(pathway)
    return map(substrates, enzyme_instances)
end

"""
    products(pathway::MetabolicPathway)

Return an `NTuple` where each entry holds the tuple of `Symbols` of product names for each enzyme, matching the
ordering produced by `enzymes(pathway)`.
"""
function products(
    pathway::MetabolicPathway,
)
    enzyme_instances = _generate_Enzymes(pathway)
    return map(products, enzyme_instances)
end

"""
    activators(pathway::MetabolicPathway)

Return an `NTuple` populated with tuple of `Symbols` of activator names for each enzyme, matching the ordering
produced by `enzymes(pathway)`.
"""
function activators(
    pathway::MetabolicPathway,
)
    enzyme_instances = _generate_Enzymes(pathway)
    return map(activators, enzyme_instances)
end

"""
    inhibitors(pathway::MetabolicPathway)

Return an `NTuple` populated with tuple of `Symbols` of inhibitor names for each enzyme, matching the ordering
produced by `enzymes(pathway)`.
"""
function inhibitors(
    pathway::MetabolicPathway,
)
    enzyme_instances = _generate_Enzymes(pathway)
    return map(inhibitors, enzyme_instances)
end

"""
    disequilibrium_ratios(pathway::MetabolicPathway, metabs::LArray, params::LArray)

Return an `NTuple` of [`disequilibrium_ratio`](@ref) values aligned with the pathway enzyme order as
reported by `enzymes(pathway)`.
"""
function disequilibrium_ratios(
    pathway::MetabolicPathway,
    metabs::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{T2,1,Vector{T2},ParamNames},
) where {T1<:Real,T2<:Real,MetabNames,ParamNames}
    enzymes = _generate_Enzymes(pathway)
    return map(enzyme -> disequilibrium_ratio(enzyme, metabs, params), enzymes)
end

"""
    reactants(pathway::MetabolicPathway)

Return the tuple of unique substrate and product names participating in the pathway reactions, ordered by first
appearance across enzymes. Excludes activators and inhibitors that are not substrates or products.
"""
@generated function reactants(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs}
    unique_reactant_names = ()
    for (e, enz) in enumerate(Enzs)
        enz_substrates = enz[2]
        enz_products = enz[3]
        for metab_name in enz_substrates
            if !(metab_name in unique_reactant_names)
                unique_reactant_names = (unique_reactant_names..., metab_name)
            end
        end
        for metab_name in enz_products
            if !(metab_name in unique_reactant_names)
                unique_reactant_names = (unique_reactant_names..., metab_name)
            end
        end
    end
    return unique_reactant_names
end

"""
    metabolites(pathway::MetabolicPathway)

Return the tuple of all unique metabolites involved in the pathway (including regulators), ordered by
first appearance.
"""
@generated function metabolites(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs}
    unique_metab_names = ()
    tuples_of_metabs = (
        substrates(MetabolicPathway{ConstMetabs,Enzs}())...,
        products(MetabolicPathway{ConstMetabs,Enzs}())...,
        activators(MetabolicPathway{ConstMetabs,Enzs}())...,
        inhibitors(MetabolicPathway{ConstMetabs,Enzs}())...,
    )
    for metab_tuple in tuples_of_metabs
        for metab_name in metab_tuple
            if !(metab_name in unique_metab_names)
                unique_metab_names = (unique_metab_names..., metab_name)
            end
        end
    end
    return unique_metab_names
end

#TODO: use labels for the enzyme and metabolite names with DimensionalData.jl
#TODO: have an option to use ConstMetabs or not
"""
    stoichiometric_matrix(pathway::MetabolicPathway)

Return the stoichiometric matrix of the metabolic pathway. Rows correspond to metabolites and columns
to enzymes ordered as in `reactants(pathway)` and `enzymes(pathway)`.
"""
@generated function stoichiometric_matrix(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs}
    metab_names = reactants(MetabolicPathway{ConstMetabs,Enzs}())
    s_matrix = zeros(Int, length(metab_names), length(Enzs))
    for (m, metab_name) in enumerate(metab_names)
        if metab_name ∉ ConstMetabs
            for (e, enz) in enumerate(Enzs)
                enz_substrates = enz[2]
                enz_products = enz[3]
                if metab_name in enz_substrates
                    stoich_coeff = sum(metab_name .== enz_substrates)
                    s_matrix[m, e] = -stoich_coeff
                end
                if metab_name in enz_products
                    stoich_coeff = sum(metab_name .== enz_products)
                    s_matrix[m, e] = stoich_coeff
                end
            end
        end
    end
    return s_matrix
end
