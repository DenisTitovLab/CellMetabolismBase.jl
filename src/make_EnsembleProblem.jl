using SciMLBase: EnsembleProblem

function make_EnsembleProblem(
    metab_path::MetabolicPathway,
    vect_init_cond::Vector{LArray},
    vect_params::Vector{LArray};
    tspan::Tuple{Float64,Float64}=(0.0, 1e8)
)
    #test that the pathway was assembled correctly
    #assert that metabs names overlap with metabolic_pathway substrates and products
    #assert that params overlap with metabolic_pathway params or maybe use the latter?

    length(vect_init_cond) == length(vect_params) || error("vect_init_cond and vect_params must have the same length")
    prob = make_ODEProblem(metab_path, vect_init_cond[1], tspan, vect_params[1])
    function prob_func(prob, i, repeat)
        @. prob.u0 = vect_init_cond[i]
        @. prob.p = vect_params[i]
        return prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, uinit=init_cond)
    return ensemble_prob
end

make_EnsembleProblem(
    metab_path::MetabolicPathway,
    vect_init_cond::Vector{LArray},
    params::LArray;
    kwargs...
) = make_EnsembleProblem(
    metab_path,
    vect_init_cond,
    [params for _ in vect_init_cond];
    kwargs...
)

make_EnsembleProblem(
    metab_path::MetabolicPathway,
    init_cond::LArray,
    vect_params::Vector{LArray};
    kwargs...
) = make_EnsembleProblem(
    metab_path,
    [init_cond for _ in vect_params],
    vect_params;
    kwargs...
)
