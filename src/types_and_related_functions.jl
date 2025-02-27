#TODO: add dependence on SciMLBase or wherever ODEProblem is defined
using LabelledArrays

#TODO: add a field to store parameter values in Enzyme

struct MetabolicPathway{ConstMetabs,Enzs} end

struct Enzyme{Name,Subs,Prods} end

@inline @generated _generate_enzymes(
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
) where {ConstMetabs,Enzs} = map(Enz -> Enzyme{Enz...}(), Enzs)

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

@inline @generated function _substrate_indxs(
    metabs::LArray{T,1,Vector{T},Syms},
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
) where {T,Syms,ConstMetabs,Enzs}
    indxs = Vector{Vector{Int}}()
    substr_names = map(Base.Fix2(getindex, 2), Enzs)
    for subsrt in substr_names
        temp_indxs = Int[]
        for m in subsrt
            if m ∉ ConstMetabs
                i = findfirst(==(m), Syms)
                push!(temp_indxs, i)
            end
        end
        push!(indxs, temp_indxs)
    end
    return indxs
end

@inline @generated function _product_indxs(
    metabs::LArray{T,1,Vector{T},Syms},
    metab_path::MetabolicPathway{ConstMetabs,Enzs},
) where {T,Syms,ConstMetabs,Enzs}
    indxs = Vector{Vector{Int}}()
    prod_names = map(Base.Fix2(getindex, 3), Enzs)
    for prod in prod_names
        temp_indxs = Int[]
        for m in prod
            if m ∉ ConstMetabs
                i = findfirst(==(m), Syms)
                push!(temp_indxs, i)
            end
        end
        push!(indxs, temp_indxs)
    end
    return indxs
end

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
