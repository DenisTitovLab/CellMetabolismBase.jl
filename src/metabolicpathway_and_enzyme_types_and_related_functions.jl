using LabelledArrays

struct MetabolicPathway{ConstMetabs,Enzs} end

struct Enzyme{Name,Subs,Prods} end

rate(enzyme::Enzyme, x, y) = error("rate function not defined for enzyme: $enzyme")

@inline @generated constant_metabs(
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = ConstMetabs

@inline @generated enzyme_names(
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Base.Fix2(getindex, 1), Enzs)

@inline @generated substrate_names(
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Base.Fix2(getindex, 2), Enzs)

@inline @generated product_names(
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Base.Fix2(getindex, 3), Enzs)

#TODO: use labels for the enzyme and metabolite names with DimensionalData.jl
#TODO: have an option to use ConstMetabs or not
@inline @generated function stoich_matrix(
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
    metabs::LArray{T,1,Vector{T},Syms},
) where {ConstMetabs,Enzs,T,Syms}
    s_matrix = zeros(Int, length(Syms), length(Enzs))
    for (i, metab_name) in enumerate(Syms)
        if metab_name ∉ ConstMetabs
            for (e, enz) in enumerate(Enzs)
                enz_substrates = enz[2]
                enz_products = enz[3]
                for m in enz_substrates
                    if m == metab_name
                        s_matrix[i, e] = -1
                    end
                end
                for m in enz_products
                    if m == metab_name
                        s_matrix[i, e] = 1
                    end
                end
            end
        end
    end
    return s_matrix
end
