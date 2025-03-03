using LabelledArrays

struct MetabolicPathway{ConstantMetabolites,Enzymes} end
function MetabolicPathway(ConstantMetabolites, Enzymes)
    ConstantMetabolites isa Tuple{Vararg{Symbol}} || error("ConstantMetabolites must be a tuple of symbols like (:Glucose, :Lactate, :ATP,)")
    Enzymes isa Tuple{Vararg{Tuple{Symbol,Tuple{Vararg{Symbol}},Tuple{Vararg{Symbol}}}}} ||
        error(
            "Second field Enzymes must be ((:Name, (:Substrate,), (:Product1,:Product2)), (:Name2, (:Substrate3,), (:Product4,)),...)",
        )
    return MetabolicPathway{ConstantMetabolites,Enzymes}()
end

struct Enzyme{Name,Substrates,Products} end
function Enzyme(Name, Substrates, Products)
    Name isa Symbol || error("Name must be a symbol like :Enz1")
    Substrates isa Tuple{Vararg{Symbol}} || error("Substrates must be a tuple of symbols like (:A_media,)")
    Products isa Tuple{Vararg{Symbol}} || error("Products must be a tuple of symbols like (:A,)")
    return Enzyme{Name,Substrates,Products}()
end

rate(enzyme::Enzyme, x, y) = error("rate function not defined for enzyme: $enzyme")

constant_metabs(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = ConstMetabs

enzyme_names(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Base.Fix2(getindex, 1), Enzs)

substrate_names(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Base.Fix2(getindex, 2), Enzs)

product_names(
    ::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Base.Fix2(getindex, 3), Enzs)

#TODO: use labels for the enzyme and metabolite names with DimensionalData.jl
#TODO: have an option to use ConstMetabs or not
@generated function stoichiometric_matrix(
    ::MetabolicPathway{ConstMetabs,Enzs},
    metabs::LArray{T,1,Vector{T},MetaboliteNames},
) where {ConstMetabs,Enzs,T,MetaboliteNames}
    s_matrix = zeros(Int, length(MetaboliteNames), length(Enzs))
    for (i, metab_name) in enumerate(MetaboliteNames)
        if metab_name ∉ ConstMetabs
            for (e, enz) in enumerate(Enzs)
                enz_substrates = enz[2]
                enz_products = enz[3]
                if metab_name in enz_substrates
                    stoich_coeff = sum(metab_name .== enz_substrates)
                    s_matrix[i, e] = -stoich_coeff
                end
                if metab_name in enz_products
                    stoich_coeff = sum(metab_name .== enz_products)
                    s_matrix[i, e] = stoich_coeff
                end
            end
        end
    end
    return s_matrix
end
