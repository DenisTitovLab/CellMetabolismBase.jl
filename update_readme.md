# CellMetabolismBase

[![Build Status](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl/graph/badge.svg?token=XC36BNU4IZ)](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

CellMetabolismBase rewrites enzyme rate equations into a system of ordinary differential equations (ODEs) that simulate the activity of metabolic pathways under various conditions in the presence and absence of allosteric regulation. The ODEs are numerically simulated using the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) library solvers in [Julia Programming Language](https://julialang.org). CellMetabolismBase package powers [CellMetabolism.jl](https://github.com/DenisTitovLab/CellMetabolism.jl) that extends CellMetabolismBase with human enzyme rate equations to enable the simulation of human cell metabolism.

## Feature Highlights

- Automatically convert enzyme `rate()` equations into `ODEProblem` and `EnsembleProblem` that can be simulated using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- `validate()` that the pathway and rate equations are correctly defined
- `remove_regulation()` of an enzyme by a regulator to facilitate the investigation of allosteric regulation in metabolic pathways
- Helper functions to calculate enzyme `rates()` and `disequilibrium_ratios()`

## Quickstart: Define -> Simulate -> Analyze

The example below tracks import, conversion, and export of a metabolite with ATP turnover. `TranspA`
brings `A` into the cell, `Enz1` combines it with `ADP` to make `B` while regenerating `ATP`,
`Enz2` converts `B` to `C`, `TranspC` exports `C`, and `ATPase` models generic ATP
consumption. `B` feeds back to inhibit `Enz1` through a Hill term with exponent two. Alongside the
regulated solution we also solve the pathway with regulation removed so you can compare both
scenarios in one figure.

**Step 1 – Describe the pathway.** We load the packages used throughout the walkthrough, define the
five-reaction network with constant media pools, and extend `rate` alongside the regulation toggle
for `Enz1`.

```julia
using CellMetabolismBase
using OrdinaryDiffEq, CairoMakie, LabelledArrays

# 1. Describe the pathway ------------------------------------------------------
pathway = MetabolicPathway(
    (:A_media, :C_media), # treated as constant
    (
        (:TranspA, (:A_media,), (:A,)),
        (:Enz1, (:A, :ADP), (:B, :ATP), (), (:B,)), # B feeds back to inhibit Enz1
        (:Enz2, (:B,), (:C,)),
        (:TranspC, (:C,), (:C_media,)),
        (:ATPase, (:ATP,), (:ADP,)),
    ),
)

function CellMetabolismBase.rate(::Enzyme{:TranspA,(:A_media,),(:A,)}, metabs, params)
    rate =
        (params.TranspA_Vmax / params.TranspA_K_A_media) *
        (metabs.A_media - metabs.A / params.TranspA_Keq) /
        (1 + metabs.A_media / params.TranspA_K_A_media + metabs.A / params.TranspA_K_A)
    return rate
end

function CellMetabolismBase.rate(
    ::Enzyme{:Enz1,(:A, :ADP),(:B, :ATP),(),(:B,)},
    metabs,
    params,
)
    rate =
        (params.Enz1_Vmax / (params.Enz1_K_A * params.Enz1_K_ADP)) *
        (metabs.A * metabs.ADP - (metabs.B * metabs.ATP) / params.Enz1_Keq) / (
            1 +
            metabs.A / params.Enz1_K_A +
            metabs.ADP / params.Enz1_K_ADP +
            metabs.B / params.Enz1_K_B +
            metabs.ATP / params.Enz1_K_ATP
        )
    inhibition_term = 1 / (1 + (metabs.B / params.Enz1_Ki_B)^2)
    return rate * inhibition_term
end

function CellMetabolismBase.remove_regulation(
    ::Enzyme{:Enz1,(:A, :ADP),(:B, :ATP),(),(:B,)},
    params,
)
    new_params = deepcopy(params)
    setproperty!(new_params, :Enz1_Ki_B, Inf)
    return new_params
end

function CellMetabolismBase.rate(::Enzyme{:Enz2,(:B,),(:C,)}, metabs, params)
    rate =
        (params.Enz2_Vmax / params.Enz2_K_B) * (metabs.B - metabs.C / params.Enz2_Keq) /
        (1 + metabs.B / params.Enz2_K_B + metabs.C / params.Enz2_K_C)
    return rate
end

function CellMetabolismBase.rate(::Enzyme{:TranspC,(:C,),(:C_media,)}, metabs, params)
    rate =
        (params.TranspC_Vmax / params.TranspC_K_C) *
        (metabs.C - metabs.C_media / params.TranspC_Keq) /
        (1 + metabs.C / params.TranspC_K_C + metabs.C_media / params.TranspC_K_C_media)
    return rate
end

function CellMetabolismBase.rate(::Enzyme{:ATPase,(:ATP,),(:ADP,)}, metabs, params)
    rate =
        (params.ATPase_Vmax / params.ATPase_K_ATP) *
        (metabs.ATP - metabs.ADP / params.ATPase_Keq) /
        (1 + metabs.ATP / params.ATPase_K_ATP + metabs.ADP / params.ATPase_K_ADP)
    return rate
end
```

**Step 2 – Configure state and parameter vectors.** We specify labelled initial conditions and the
parameter set used by the rate laws.

```julia
# 2. Initial conditions and parameters -----------------------------------------
metabs0 = LVector(
    A_media = 5.0,
    C_media = 0.0,
    A = 0.1,
    B = 0.0,
    C = 0.0,
    ATP = 1.5,
    ADP = 1.0,
)
params = LVector(
    TranspA_Vmax = 1.0,
    TranspA_Keq = 8.0,
    TranspA_K_A_media = 1.0,
    TranspA_K_A = 1.0,
    Enz1_Vmax = 0.8,
    Enz1_Keq = 15.0,
    Enz1_K_A = 0.4,
    Enz1_K_ADP = 0.3,
    Enz1_K_B = 0.3,
    Enz1_K_ATP = 0.5,
    Enz1_Ki_B = 0.6,
    Enz2_Vmax = 1.2,
    Enz2_Keq = 6.0,
    Enz2_K_B = 0.4,
    Enz2_K_C = 0.5,
    TranspC_Vmax = 1.0,
    TranspC_Keq = 5.0,
    TranspC_K_C = 0.5,
    TranspC_K_C_media = 0.5,
    ATPase_Vmax = 0.6,
    ATPase_Keq = 20.0,
    ATPase_K_ATP = 0.7,
    ATPase_K_ADP = 0.6,
)
```

**Step 3 – Build problems and solve both scenarios.** We create the regulated ODE problem, use the
`remove_regulation` helper for `Enz1`, and solve a second problem with feedback disabled.

```julia
# 3. Build and solve the ODE problems ------------------------------------------
tspan = (0.0, 40.0)
prob = make_ODEProblem(pathway, metabs0, tspan, params)
sol = solve(prob, RadauIIA9(); saveat = range(tspan...; length = 300), abstol = 1e-15, reltol = 1e-8)

enz1 = CellMetabolismBase.Enzyme(:Enz1, (:A, :ADP), (:B, :ATP), (), (:B,))
params_noreg = remove_regulation(enz1, params)
prob_noreg = make_ODEProblem(pathway, metabs0, tspan, params_noreg)
sol_noreg = solve(prob_noreg, RadauIIA9(); saveat = sol.t, abstol = 1e-15, reltol = 1e-8)
```

**Step 4 – Compute diagnostics and visualise the outcomes.** We compare metabolites, enzyme fluxes,
and disequilibrium ratios for both trajectories. Dotted lines reuse the solid colours so you can see
the impact of removing regulation at a glance.

```julia
# 4. Analyse metabolite concs, rates, and disequilibrium ratios ---------------
concs = sol.u
concs_noreg = sol_noreg.u
rates_over_time = [rates(pathway, state, params) for state in concs]
rates_over_time_noreg = [rates(pathway, state, params_noreg) for state in concs_noreg]
ratio_over_time = [disequilibrium_ratios(pathway, state, params) for state in concs]
ratio_over_time_noreg =
    [disequilibrium_ratios(pathway, state, params_noreg) for state in concs_noreg]
enzyme_labels = enzymes(pathway)

fig = Figure(size = (960, 320))
colors = Makie.wong_colors()

ax_concs = Axis(fig[1, 1], title = "Metabolites", xlabel = "time", ylabel = "concentration")
for (i, metab) in enumerate(propertynames(metabs0))
    color = colors[mod1(i, length(colors))]
    lines!(
        ax_concs,
        sol.t,
        [s[metab] for s in concs];
        color = color,
        label = string(metab),
    )
    lines!(
        ax_concs,
        sol_noreg.t,
        [s[metab] for s in concs_noreg];
        color = color,
        linestyle = :dot,
    )
end
axislegend(ax_concs, position = :rb)

ax_rates = Axis(fig[1, 2], title = "Rates()", xlabel = "time", ylabel = "flux")
for (j, label) in enumerate(enzyme_labels)
    color = colors[mod1(j, length(colors))]
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
axislegend(ax_rates, position = :rb)

ax_ratio = Axis(fig[1, 3], title = "Disequilibrium", xlabel = "time", ylabel = "Q / K_eq")
for (j, label) in enumerate(enzyme_labels)
    color = colors[mod1(j, length(colors))]
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
axislegend(ax_ratio, position = :rb)

display(fig)
```

The solid traces correspond to the regulated pathway, while dotted lines show the same colour for
the deregulated counterpart. Removing `B`’s feedback on `Enz1` accelerates `A` consumption, raises
`B` and `C` production, and shifts fluxes in a way that keeps disequilibrium ratios farther from
unity—an effective comparison for reports or teaching material.

## Roadmap

- Add Global Sensitivity Analysis functionality to identify key parameters influencing pathway behavior
- Introduce unit-aware parameters and initial conditions to ensure unit consistency
- Implement atom conservation validation to catch errors in metabolic reaction formulation
- Implement ability to write down isotope tracing ODE models

## CellMetabolismBase vs Catalyst.jl

[Catalyst.jl](https://catalyst.sciml.ai/stable/) is a mature package with many features for analysis and simulation of chemical reaction networks in [Julia Programming Language](https://julialang.org) and you should use it instead of CellMetabolismBase for most applications. The specific focus of CellMetabolismBase is to facilitate the investigation of control of metabolic pathways by allosteric regulators while [Catalyst.jl](https://catalyst.sciml.ai/stable/) is a more general framework for chemical reaction network modeling. The goal behind CellMetabolismBase is to have a lightweight package with small number of dependencies that converts a collection of enzyme rate equations into a metabolic pathway that can be simulated with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) and provides utilities that validate that enzyme rate equations and metabolic pathways are written down correctly.  CellMetabolismBase is designed to be extended by other packages such as [CellMetabolism.jl](https://github.com/DenisTitovLab/CellMetabolism.jl) to enable simulation of large cellular metabolism models with custom rate equations and regulation patterns.
