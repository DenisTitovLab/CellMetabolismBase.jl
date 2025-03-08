using SciMLBase: EnsembleProblem
using Distributions

"""
    make_EnsembleProblem(
    metab_path::MetabolicPathway,
    vect_init_cond::Vector{LArray{T1,1,Vector{T1},MetabNames}},
    vect_params::Vector{LArray{T2,1,Vector{T2},ParamNames}};
    tspan::Tuple{T3,T3}=(0.0, 1e8)
) where {T1<:Real,T2<:Real,T3<:Real,MetabNames,ParamNames}

Create an `EnsembleProblem` for a metabolic pathway simulation with multiple initial conditions and parameters.

# Arguments
- `metab_path::MetabolicPathway`: The metabolic pathway model to simulate.
- `vect_init_cond::Vector{<:LArray}`: Vector of initial conditions for each ensemble member.
- `vect_params::Vector{<:LArray}`: Vector of parameters for each ensemble member.
- `tspan=(0.0, 1e8)`: Time span for the simulation.

# Returns
- `EnsembleProblem`: An ensemble problem that can be solved using DifferentialEquations.jl's ensemble solvers.

# Notes
- `vect_init_cond` and `vect_params` must have the same length.
- Each element in the ensemble will use the corresponding initial conditions and parameters from the provided vectors.
- The ensemble problem is constructed from a base ODE problem created with the first elements of the initial conditions and parameters vectors.

# Example
```julia
using CellMetabolismBase
using DifferentialEquations

metab_path = MetabolicPathway(...)
init_cond = LArray(...)
params = LArray(...)
tspan = (0.0, 1e8)
ensemble_prob = make_EnsembleProblem(metab_path, [init_cond1, init_cond2], [params1, params2], tspan)
ensemble_sol = solve(ensemble_prob, RadauIIA9(), trajectories=2)
```
"""
function make_EnsembleProblem(
    metab_path::MetabolicPathway,
    vect_init_cond::Vector{LArray{T1,1,Vector{T1},MetabNames}},
    vect_params::Vector{LArray{T2,1,Vector{T2},ParamNames}};
    tspan::Tuple{T3,T3}=(0.0, 1e8)
) where {T1<:Real,T2<:Real,T3<:Real,MetabNames,ParamNames}
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
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    return ensemble_prob
end

make_EnsembleProblem(
    metab_path::MetabolicPathway,
    vect_init_cond::Vector{LArray{T1,1,Vector{T1},MetabNames}},
    params::LArray{T2,1,Vector{T2},ParamNames};
    kwargs...
) where {T1<:Real,T2<:Real,MetabNames,ParamNames} = make_EnsembleProblem(
    metab_path,
    vect_init_cond,
    [params for _ in vect_init_cond];
    kwargs...
)

make_EnsembleProblem(
    metab_path::MetabolicPathway,
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    vect_params::Vector{LArray{T2,1,Vector{T2},ParamNames}};
    kwargs...
) where {T1<:Real,T2<:Real,MetabNames,ParamNames} = make_EnsembleProblem(
    metab_path,
    [init_cond for _ in vect_params],
    vect_params;
    kwargs...
)

make_EnsembleProblem(
    metab_path::MetabolicPathway,
    init_cond::LArray{<:Distribution,1,<:Vector{<:Distribution},MetabNames},
    params::LArray{<:Distribution,1,<:Vector{<:Distribution},ParamNames};
    n_bootstraps::Int=1000,
    kwargs...
) where {MetabNames,ParamNames} = make_EnsembleProblem(
    metab_path,
    [rand.(init_cond) for _ in 1:n_bootstraps],
    [rand.(params) for _ in 1:n_bootstraps];
    kwargs...
)


make_EnsembleProblem(
    metab_path::MetabolicPathway,
    init_cond::LArray{<:Distribution,1,<:Vector{<:Distribution},MetabNames},
    params::LArray{T,1,Vector{T},ParamNames};
    n_bootstraps::Int=1000,
    kwargs...
) where {T<:Real,MetabNames,ParamNames} = make_EnsembleProblem(
    metab_path,
    [rand.(init_cond) for _ in 1:n_bootstraps],
    params;
    kwargs...
)

make_EnsembleProblem(
    metab_path::MetabolicPathway,
    init_cond::LArray{T1,1,Vector{T1},MetabNames},
    params::LArray{<:Distribution,1,<:Vector{<:Distribution},ParamNames};
    n_bootstraps::Int=1000,
    kwargs...
) where {T1<:Real,MetabNames,ParamNames} = make_EnsembleProblem(
    metab_path,
    init_cond,
    [rand.(params) for _ in 1:n_bootstraps];
    kwargs...
)
