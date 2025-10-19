# CellMetabolismBase

[![Build Status](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl/graph/badge.svg?token=XC36BNU4IZ)](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl)

CellMetabolismBase turns symbolic pathway descriptions into DifferentialEquations.jl models so you can explore steady states, dynamics, and regulation patterns with minimal boilerplate. The package powers CellMetabolism.jl and is designed to be extended with custom rate laws, regulation toggles, and analysis tooling.

## Quickstart: Define -> Simulate -> Analyze

The example below walks through a three-enzyme pathway with feedback inhibition (metabolite `B` inhibits `Enz1`). We simulate the dynamics, compute instantaneous rates and disequilibrium ratios, and visualize everything side-by-side using CairoMakie.

```julia
using CellMetabolismBase
using OrdinaryDiffEq, CairoMakie, LabelledArrays

# 1. Describe the pathway ------------------------------------------------------
pathway = MetabolicPathway(
    (:A_media,), # treated as constant
    (
        (:Enz1, (:A_media,), (:A,), (), (:B,)), # B feeds back to inhibit Enz1
        (:Enz2, (:A,), (:B,)),
        (:Enz3, (:B,), (:C,)),
    ),
)

# Extend CellMetabolismBase.rate for each enzyme. Real projects usually place
# these in src/ files, but they work the same in scripts or Pluto notebooks.
function CellMetabolismBase.rate(
    ::Enzyme{:Enz1,(:A_media,),(:A,),(),(:B,)},
    metabs,
    params,
)
    forcing = params.Enz1_Vmax * (metabs.A_media - metabs.A / params.Enz1_Keq)
    inhibition = 1 / (1 + metabs.B / params.Enz1_Ki_B)
    return forcing * inhibition
end

function CellMetabolismBase.rate(::Enzyme{:Enz2,(:A,),(:B,)}, metabs, params)
    drive = params.Enz2_Vmax * (metabs.A - metabs.B / params.Enz2_Keq)
    sat = 1 / (1 + metabs.A / params.Enz2_K_A + metabs.B / params.Enz2_K_B)
    return drive * sat
end

function CellMetabolismBase.rate(::Enzyme{:Enz3,(:B,),(:C,)}, metabs, params)
    drive = params.Enz3_Vmax * (metabs.B - metabs.C / params.Enz3_Keq)
    return drive / (1 + metabs.B / params.Enz3_K_B)
end

# 2. Initial conditions and parameters -----------------------------------------
metabs0 = LVector(A_media = 5.0, A = 0.0, B = 0.0, C = 0.0)
params = LVector(
    Enz1_Vmax = 1.0,
    Enz1_Keq = 8.0,
    Enz1_Ki_B = 0.6,
    Enz2_Vmax = 1.4,
    Enz2_Keq = 3.0,
    Enz2_K_A = 0.4,
    Enz2_K_B = 0.7,
    Enz3_Vmax = 2.0,
    Enz3_Keq = 4.0,
    Enz3_K_B = 0.5,
)

# 3. Build and solve the ODE problem -------------------------------------------
tspan = (0.0, 40.0)
prob = make_ODEProblem(pathway, metabs0, tspan, params)
sol = solve(prob, Tsit5(); saveat = range(tspan...; length = 300))

# 4. Analyse metabolite states, rates, and disequilibrium ratios ---------------
states = sol.u
rates_over_time = [rates(pathway, state, params) for state in states]
ratio_over_time = [disequilibrium_ratios(pathway, state, params) for state in states]
enzyme_labels = enzymes(pathway)

fig = Figure(resolution = (960, 320))
colors = Makie.wong_colors()

ax_states = Axis(fig[1, 1], title = "Metabolites", xlabel = "time", ylabel = "concentration")
for (i, metab) in enumerate(propertynames(metabs0))
    lines!(ax_states, sol.t, [s[metab] for s in states], color = colors[i], label = string(metab))
end
axislegend(ax_states, position = :rb)

ax_rates = Axis(fig[1, 2], title = "Rates()", xlabel = "time", ylabel = "flux")
for (j, label) in enumerate(enzyme_labels)
    lines!(ax_rates, sol.t, getindex.(rates_over_time, j), color = colors[j], label = String(label))
end

ax_ratio = Axis(fig[1, 3], title = "Disequilibrium", xlabel = "time", ylabel = "Q / K_eq")
for (j, label) in enumerate(enzyme_labels)
    lines!(ax_ratio, sol.t, getindex.(ratio_over_time, j), color = colors[j], label = String(label))
    hlines!(ax_ratio, [1.0]; color = :gray80, linestyle = :dash)
end

fig
```

This compact workflow highlights how CellMetabolismBase keeps pathway definitions, integrators, and diagnostics in sync. The three panels share the same color coding so you can connect metabolite behaviour to flux changes and equilibrium shifts at a glance.

### Remove Feedback and Compare

Regulation toggles integrate naturally with the `remove_regulation` hooks that CellMetabolismBase expects downstream packages to extend. For the pathway above we can neutralize the feedback by redefining the inhibitor constant:

```julia
function CellMetabolismBase.remove_regulation(
    ::Enzyme{:Enz1,(:A_media,),(:A,),(),(:B,)},
    params,
)
    new_params = deepcopy(params)
    setproperty!(new_params, :Enz1_Ki_B, Inf)
    return new_params
end

params_noreg = remove_regulation(first(enzymes(pathway)), params)
sol_noreg = solve(make_ODEProblem(pathway, metabs0, tspan, params_noreg), Tsit5(); saveat = sol.t)
```

Plotting `sol` and `sol_noreg` together reveals how removing feedback accelerates `A` consumption and shifts disequilibrium ratios beyond unity - an effective demonstration for reports or teaching material.

## Feature Highlights

| Task | Helper | Notes |
| --- | --- | --- |
| Build ODE problems | `make_ODEProblem(pathway, metabs0, tspan, params)` | Returns fully typed SciML problems. |
| Explore variability | `make_EnsembleProblem(pathway, initial_conditions, params; tspan)` | Works with any DifferentialEquations.jl ensemble algorithm. |
| Inspect instantaneous flux | `rates(pathway, state, params)` | Ordered to match `enzymes(pathway)`; great for plotting or control analysis. |
| Measure distance from equilibrium | `disequilibrium_ratios(pathway, state, params)` | Equals one at thermodynamic equilibrium. |
| Toggle regulation | Extend `remove_regulation` | Override per regulator or all-at-once to study knockouts. |

## Beyond the Quickstart

- **Pathway variants:** keep multiple `MetabolicPathway` definitions side-by-side to compare regulation topologies or substrate choices.
- **Sensitivity analysis:** plug the generated problems into Global Sensitivity Analysis (GSA) tooling from SciML to discover leverage points.
- **Performance tips:** rates are generic over `LArray` inputs; use `StaticArrays`-backed structures or warm-started integrators for tight loops.
- **Documentation build:** `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include(\"docs/make.jl\")'` regenerates the docs site (Documenter.jl).

## Contributing Checklist

- Instantiate dependencies once per clone: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- Run tests (Aqua, JET, and feature suites): `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Format touched files with `JuliaFormatter.format(".", indent = 4)` before submitting PRs.
- Summarize functional changes and call out API updates in your PR description - link to relevant figures or notebooks when possible.

## Additional Resources

- The `test/` directory bundles minimal pathways and helper fixtures; copy them for bespoke experiments.
- `docs/src/tutorials/` covers more advanced workflows such as isotope tracing and batch-mode ensembles.
- For real-world enzyme kinetics, see [CellMetabolism.jl](https://github.com/DenisTitovLab/CellMetabolism.jl), which extends this base layer with curated rate laws.
