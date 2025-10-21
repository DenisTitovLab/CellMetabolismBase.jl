# CellMetabolismBase

[![Build Status](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl/graph/badge.svg?token=XC36BNU4IZ)](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

CellMetabolismBase rewrites enzyme rate equations into ordinary differential equation (ODE) models that can be used to simulate the activity of metabolic pathways under various conditions in the presence and absence of allosteric regulation. The ODE models are numerically simulated using the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) library solvers in [Julia Programming Language](https://julialang.org). CellMetabolismBase package powers [CellMetabolism.jl](https://github.com/DenisTitovLab/CellMetabolism.jl) that extends CellMetabolismBase with human enzyme rate equations to enable the simulation of human cell metabolism.

## Feature Highlights

- Automatically convert enzyme `rate()` equations into `ODEProblem` and `EnsembleProblem` that can be simulated using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- `validate()` that the pathway and rate equations are correctly defined
- `remove_regulation()` of an enzyme by a regulator to facilitate the investigation of allosteric regulation in metabolic pathways
- Helper functions to calculate enzyme `rates()` and `disequilibrium_ratios()`

## Quickstart: Define -> Simulate -> Analyze

The example below revisits the three-enzyme ABC pathway from the README. `Enz2` converts `A` into
two copies of `B` and is inhibited by `B` through a Hill term with exponent two. Alongside the
regulated solution we also solve the pathway with regulation removed so you can compare both
scenarios in one figure.

**Step 1 – Describe the pathway.** We load the packages used throughout the walkthrough, define the
ABC network, and extend `rate` with kinetics that include Hill-type feedback on `Enz2`.

```julia
using CellMetabolismBase
using OrdinaryDiffEq, CairoMakie, LabelledArrays

# 1. Describe the pathway ------------------------------------------------------
pathway = MetabolicPathway(
    (:A_media,), # treated as constant
    (
        (:Enz1, (:A_media,), (:A,)),
        (:Enz2, (:A,), (:B, :B), (), (:B,)), # B feeds back to inhibit Enz2
        (:Enz3, (:B,), (:C,)),
    ),
)

function CellMetabolismBase.rate(::Enzyme{:Enz1,(:A_media,),(:A,)}, metabs, params)
    return params.Enz1_Vmax *
           (metabs.A_media - metabs.A / params.Enz1_Keq) /
           (1 + metabs.A_media / params.Enz1_K_A_media + metabs.A / params.Enz1_K_A)
end

function CellMetabolismBase.rate(::Enzyme{:Enz2,(:A,),(:B, :B),(),(:B,)}, metabs, params)
    drive =
        params.Enz2_Vmax *
        (metabs.A - metabs.B^2 / params.Enz2_Keq) /
        (1 + metabs.A / params.Enz2_K_A + (metabs.B / params.Enz2_K_B)^2)
    inhibition = 1 / (1 + (metabs.B / params.Enz2_Ki_B)^2)
    return drive * inhibition
end

function CellMetabolismBase.remove_regulation(
    ::Enzyme{:Enz2,(:A,),(:B, :B),(),(:B,)},
    params,
)
    new_params = deepcopy(params)
    setproperty!(new_params, :Enz2_Ki_B, Inf)
    return new_params
end

function CellMetabolismBase.rate(::Enzyme{:Enz3,(:B,),(:C,)}, metabs, params)
    return params.Enz3_Vmax *
           (metabs.B - metabs.C / params.Enz3_Keq) /
           (1 + metabs.B / params.Enz3_K_B + metabs.C / params.Enz3_Keq)
end
```

**Step 2 – Configure state and parameter vectors.** We specify labelled initial conditions and the
parameter set used by the rate laws.

```julia
# 2. Initial conditions and parameters -----------------------------------------
metabs0 = LVector(A_media = 5.0, A = 0.0, B = 0.0, C = 0.0)
params = LVector(
    Enz1_Vmax = 1.0,
    Enz1_Keq = 10.0,
    Enz1_K_A_media = 1.0,
    Enz1_K_A = 1.0,
    Enz2_Vmax = 1.5,
    Enz2_Keq = 5.0,
    Enz2_K_A = 0.5,
    Enz2_K_B = 0.5,
    Enz2_Ki_B = 0.6,
    Enz3_Vmax = 2.0,
    Enz3_Keq = 4.0,
    Enz3_K_B = 0.3,
)
```

**Step 3 – Build problems and solve both scenarios.** We create the regulated ODE problem, extend
`remove_regulation` for `Enz2`, and solve a second problem with feedback disabled.

```julia
# 3. Build and solve the ODE problems ------------------------------------------
tspan = (0.0, 40.0)
prob = make_ODEProblem(pathway, metabs0, tspan, params)
sol = solve(prob, RadauIIA9(); saveat = range(tspan...; length = 300), abstol = 1e-15, reltol = 1e-8)

enz2 = CellMetabolismBase.Enzyme(:Enz2, (:A,), (:B, :B), (), (:B,))
params_noreg = remove_regulation(enz2, params)
prob_noreg = make_ODEProblem(pathway, metabs0, tspan, params_noreg)
sol_noreg = solve(prob_noreg, RadauIIA9(); saveat = sol.t, abstol = 1e-15, reltol = 1e-8)
```

**Step 4 – Compute diagnostics and visualise the outcomes.** We compare metabolites, enzyme fluxes,
and disequilibrium ratios for both trajectories. Dotted lines reuse the solid colours so you can see
the impact of removing regulation at a glance.

```julia
# 4. Analyse metabolite states, rates, and disequilibrium ratios ---------------
states = sol.u
states_noreg = sol_noreg.u
rates_over_time = [rates(pathway, state, params) for state in states]
rates_over_time_noreg = [rates(pathway, state, params_noreg) for state in states_noreg]
ratio_over_time = [disequilibrium_ratios(pathway, state, params) for state in states]
ratio_over_time_noreg =
    [disequilibrium_ratios(pathway, state, params_noreg) for state in states_noreg]
enzyme_labels = enzymes(pathway)

fig = Figure(size = (960, 320))
colors = Makie.wong_colors()

ax_states = Axis(fig[1, 1], title = "Metabolites", xlabel = "time", ylabel = "concentration")
for (i, metab) in enumerate(propertynames(metabs0))
    color = colors[i]
    lines!(
        ax_states,
        sol.t,
        [s[metab] for s in states];
        color = color,
        label = string(metab),
    )
    lines!(
        ax_states,
        sol_noreg.t,
        [s[metab] for s in states_noreg];
        color = color,
        linestyle = :dot,
    )
end
axislegend(ax_states, position = :rb)

ax_rates = Axis(fig[1, 2], title = "Rates()", xlabel = "time", ylabel = "flux")
for (j, label) in enumerate(enzyme_labels)
    color = colors[j]
    lines!(
        ax_rates,
        sol.t,
        getindex.(rates_over_time, j);
        color = color,
        label = String(label),
    )
    lines!(
        ax_rates,
        sol_noreg.t,
        getindex.(rates_over_time_noreg, j);
        color = color,
        linestyle = :dot,
    )
end

ax_ratio = Axis(fig[1, 3], title = "Disequilibrium", xlabel = "time", ylabel = "Q / K_eq")
for (j, label) in enumerate(enzyme_labels)
    color = colors[j]
    lines!(
        ax_ratio,
        sol.t,
        getindex.(ratio_over_time, j);
        color = color,
        label = String(label),
    )
    lines!(
        ax_ratio,
        sol_noreg.t,
        getindex.(ratio_over_time_noreg, j);
        color = color,
        linestyle = :dot,
    )
    hlines!(ax_ratio, [1.0]; color = :gray80, linestyle = :dash)
end

fig
```

The solid traces correspond to the regulated pathway, while dotted lines show the same colour for
the deregulated counterpart. Removing `B`’s feedback accelerates `A` consumption, raises `B`
concentrations, and shifts fluxes in a way that keeps disequilibrium ratios farther from unity—an
effective comparison for reports or teaching material.

## Roadmap

- Introduce unit-aware parameters and initial conditions to ensure dimensional consistency
- Add Global Sensitivity Analysis functionality to identify key parameters influencing pathway behavior
- Implement ability to write down isotope tracing ODE models
- Implement atom conservation validation to catch errors in metabolic pathway definitions

## CellMetabolismBase vs Catalyst.jl

[Catalyst.jl](https://catalyst.sciml.ai/stable/) is the mature package with many features for analysis and simulation of chemical reaction networks in [Julia Programming Language](https://julialang.org) and you should use it instead of CellMetabolismBase for most applications. The specific focus of CellMetabolismBase is to facilitate the investigation of control of metabolic pathways by allosteric regulators while [Catalyst.jl](https://catalyst.sciml.ai/stable/) is a more general framework for chemical reaction network modeling. The goal behind CellMetabolismBase is to have a lightweight package with small number of dependencies that converts a collection of enzyme rate equations into a metabolic pathway that can be simulated with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) and provides utilities that ensure that many typos in enzyme rate equations and pathway formulation are caught early.  CellMetabolismBase is designed to be extended by other packages such as [CellMetabolism.jl](https://github.com/DenisTitovLab/CellMetabolism.jl) to enable simulation of large cellular metabolism models with custom rate equations and regulation patterns.
