using LabelledArrays

"""
    MetabolicPathway{ConstantMetabolites,Enzymes}

Type-level description of a metabolic pathway where metabolite and enzyme names are embedded in the
type parameters. Construct instances via [`MetabolicPathway`](@ref).
"""
struct MetabolicPathway{ConstantMetabolites,Enzymes} end

#TODO: transition to using DSL and move all the synthax checking to the DSL
"""
    MetabolicPathway(ConstantMetabolites, Enzymes)

Validate the provided constant metabolite symbols and enzyme tuples before constructing a
`MetabolicPathway`. `ConstantMetabolites` must be a tuple of symbols and `Enzymes` a tuple of tuples in
the format accepted by [`Enzyme`](@ref).
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

Type-level representation of an enzyme with named substrates, products, activators, and inhibitors.
Construct instances via [`Enzyme`](@ref).
"""
struct Enzyme{Name,Substrates,Products,Activators,Inhibitors} end

#TODO: transition to using DSL and move all the synthax checking to the DSL
"""
    Enzyme(Name, Substrates, Products[, Activators, Inhibitors])

Validate metabolite and regulator tuples before constructing an `Enzyme`. `Name` must be a symbol,
`Substrates` and `Products` are tuples of symbols with at least one entry, and `Activators` and
`Inhibitors` are tuples of symbols (possibly empty).
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
    enzyme_rate(enzyme::Enzyme, metabs, params)

Fallback method for enzyme rate calculations. Packages depending on `CellMetabolismBase` should extend
this function for relevant enzyme variants. The default implementation throws an error to highlight
missing rate equations.
"""
enzyme_rate(enzyme::Enzyme, metabs, params) = error("rate function not defined for enzyme: $enzyme")

"""
    enzyme_name(enzyme::Enzyme)

Return `Symbol` of enzyme name.
"""
@inline enzyme_name(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Name

"""
    substrates_name(enzyme::Enzyme)

Return the tuple of `Symbols` of substrate names of the enzyme reaction.
"""
@inline substrates_name(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Substrates

"""
    products_name(enzyme::Enzyme)

Return the tuple of `Symbols` of product names of the enzyme reaction.
"""
@inline products_name(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Products

"""
    activators_name(enzyme::Enzyme)

Return the tuple of `Symbols` of activator names of the enzyme.
"""
@inline activators_name(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Activators

"""
    inhibitors_name(enzyme::Enzyme)

Return the tuple of `Symbols` of inhibitor names of the enzyme.
"""
@inline inhibitors_name(
    ::Enzyme{Name,Substrates,Products,Activators,Inhibitors},
) where {Name,Substrates,Products,Activators,Inhibitors} = Inhibitors

"""
    disequilibrium_ratio(enzyme::Enzyme, metabs::LArray, params::LArray)

Compute the ratio of product of product concentrations to product of substrate concentrations scaled by the equilibrium constant (i.e., disequilibrium ratio) for the enzyme. Throws an `ArgumentError` if the corresponding `_Keq` parameter is missing in `params`.
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
    constant_metabs(pathway::MetabolicPathway)

Return the tuple of `Symbol` of metabolite names treated as constants in the pathway definition.
"""
constant_metabs(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = ConstMetabs

"""
    enzyme_names(pathway::MetabolicPathway)

Return an `NTuple` with `Symbols` of enzyme names in pathway order.
"""
function enzyme_names(
    pathway::MetabolicPathway,
)
    enzymes = _generate_Enzymes(pathway)
    return map(enzyme_name, enzymes)
end

"""
    substrates_names(pathway::MetabolicPathway)

Return an `NTuple` where each entry holds the tuple of `Symbols` of substrate names for each enzyme, matching
the ordering produced by `enzyme_names(pathway)`.
"""
function substrates_names(
    pathway::MetabolicPathway,
)
    enzymes = _generate_Enzymes(pathway)
    return map(substrates_name, enzymes)
end

"""
    products_names(pathway::MetabolicPathway)

Return an `NTuple` where each entry holds the tuple of `Symbols` of product names for each enzyme, matching the
ordering produced by `enzyme_names(pathway)`.
"""
function products_names(
    pathway::MetabolicPathway,
)
    enzymes = _generate_Enzymes(pathway)
    return map(products_name, enzymes)
end

"""
    activators_names(pathway::MetabolicPathway)

Return an `NTuple` populated with tuple of `Symbols` of activator names for each enzyme, matching the ordering
produced by `enzyme_names(pathway)`.
"""
function activators_names(
    pathway::MetabolicPathway,
)
    enzymes = _generate_Enzymes(pathway)
    return map(activators_name, enzymes)
end

"""
    inhibitors_names(pathway::MetabolicPathway)

Return an `NTuple` populated with tuple of `Symbols` of inhibitor names for each enzyme, matching the ordering
produced by `enzyme_names(pathway)`.
"""
function inhibitors_names(
    pathway::MetabolicPathway,
)
    enzymes = _generate_Enzymes(pathway)
    return map(inhibitors_name, enzymes)
end

"""
    disequilibrium_ratios(pathway::MetabolicPathway, metabs::LArray, params::LArray)

Return an `NTuple` of [`disequilibrium_ratio`](@ref) values aligned with the pathway enzyme order as
reported by `enzyme_names(pathway)`.
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
    reactant_names(pathway::MetabolicPathway)

Return the tuple of unique substrate and product names participating in the pathway reactions, ordered by first
appearance across enzymes. Excludes activators and inhibitors that are not substrates or products.
"""
@generated function reactant_names(
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
    all_metabolite_names(pathway::MetabolicPathway)

Return the tuple of all unique metabolites involved in the pathway (including regulators), ordered by
first appearance.
"""
@generated function all_metabolite_names(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs}
    unique_metab_names = ()
    tuples_of_metabs = (
        substrates_names(MetabolicPathway{ConstMetabs,Enzs}())...,
        products_names(MetabolicPathway{ConstMetabs,Enzs}())...,
        activators_names(MetabolicPathway{ConstMetabs,Enzs}())...,
        inhibitors_names(MetabolicPathway{ConstMetabs,Enzs}())...,
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
to enzymes ordered as in `reactant_names(pathway)` and `enzyme_names(pathway)`.
"""
@generated function stoichiometric_matrix(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs}
    metab_names = reactant_names(MetabolicPathway{ConstMetabs,Enzs}())
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
