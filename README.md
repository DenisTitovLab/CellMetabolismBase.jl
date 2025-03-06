# CellMetabolismBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DenisTitovLab.github.io/CellMetabolismBase.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DenisTitovLab.github.io/CellMetabolismBase.jl/dev/)
[![Build Status](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/CellMetabolismBase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl/graph/badge.svg?token=XC36BNU4IZ)](https://codecov.io/gh/DenisTitovLab/CellMetabolismBase.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

CellMetabolismBase.jl is a framework for simulating and analyzing cellular metabolic pathways. The goal of the package is to provide a convenient interface to automatically convert a list of enzymes into a form that can be used by ordinary differential equations (ODE) solvers of [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) and [Scientific Machine Learning (SciML)](https://sciml.ai) ecosystems in [Julia Programming Language](https://julialang.org). CellMetabolismBase.jl powers CellMetabolism.jl that enables simulation and analysis of human cell metabolism using a experimentally determined enzyme rate equations.

## Features

- Define `MetabolicPathway`s using enzyme rate equations
- Automatically create DifferentialEquations.jl-compatible ODE models from `MetabolicPathway` definitions
- Simulate `MetabolicPathway` activity using DifferentialEquations.jl at a wide range of parameter values and initial conditions
- Estimate model prediction uncertainty for comparison with experimental data
- Automated discovery of parameters that control a specific pathway behaviour using Global Sensitivity Analysis (GSA)

## Roadmap

- Add functionality to validate that `MetabolicPathway` is defined correctly, such as testing enzyme rate equations, determining conserved moieties, and ensuring conservation of atoms.
- Add functionality for Global Sensitivity Analysis (GSA) to identify parameters that control a specific pathway behaviour.

## Installation

```julia
using Pkg
Pkg.add("CellMetabolismBase")
```

Or in pkg mode (press `]` in REPL):

```
add CellMetabolismBase
```

## Basic Example

```julia
using CellMetabolismBase, OrdinaryDiffEq, LabelledArrays, CairoMakie

# We'll investigate a simple metabolic pathway with three enzymes:
# Enz1: A_media → A
# Enz2: A → 2B
# Enz3: B → C

# Define a simple enzyme-catalyzed pathway consisting of three enzymes:
# First argument is a tuple of constant metabolites (if any)
# Second argument is a tuple of enzyme definitions
# Each enzyme is defined as (enzyme_name, (substrates...), (products...))
pathway = MetabolicPathway(
    (:A_media,),
    ((:Enz1, (:A_media,), (:A,)), (:Enz2, (:A,), (:B, :B)), (:Enz3, (:B,), (:C,))),
)

# Define enzyme rate laws
function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz1,(:A_media,),(:A,)}, metabs, params)
    rate =
        params.Enz1_Vmax * (metabs.A_media - metabs.A / params.Enz1_Keq) /
        (1 + metabs.A_media / params.Enz1_K_A_media + metabs.A / params.Enz1_K_A)
    return rate
end
function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz2,(:A,),(:B, :B)}, metabs, params)
    rate =
        params.Enz2_Vmax * (metabs.A - metabs.B^2 / params.Enz2_Keq) /
        (1 + metabs.A / params.Enz2_K_A + (metabs.B / params.Enz2_K_B)^2)
    return rate
end
function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz3,(:B,),(:C,)}, metabs, params)
    rate =
        params.Enz3_Vmax * (metabs.B - metabs.C / params.Enz3_Keq) /
        (1 + metabs.B / params.Enz3_K_B + metabs.C / params.Enz3_Keq)
    return rate
end

# Set initial conditions
init_cond = LVector(A_media = 5.0, A = 0.0, B = 0.0, C = 0.0)

# Set parameters
params = LVector(
    Enz1_Vmax = 1.0,
    Enz1_Keq = 10.0,
    Enz1_K_A_media = 1.0,
    Enz1_K_A = 1.0,
    Enz2_Vmax = 1.5,
    Enz2_Keq = 5.0,
    Enz2_K_A = 0.5,
    Enz2_K_B = 0.5,
    Enz3_Vmax = 2.0,
    Enz3_Keq = 4.0,
    Enz3_K_B = 0.3,
)

# Create and solve ODE problem
tspan = (0.0, 10000.0)
prob = make_ODEProblem(pathway, init_cond, tspan, params)
sol = solve(prob, Rodas5P(), abstol = 1e-15, reltol = 1e-8)

# Plot results
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Concentration")
for metabolite in propertynames(init_cond)
    lines!(ax, sol.t, [u[metabolite] for u in sol.u], label = string(metabolite))
end
axislegend(ax)
fig
```
## Ensemble Simulations

Ensemble simulations allow you to explore parameter or initial condition variations:

```julia
using CellMetabolismBase
using DifferentialEquations
using LabelledArrays

# Define initial condition and parameters
init_cond = LVector(A_media=2.0, A=0.5, B=0.5, C=0.5, D=0.0)

# Define multiple parameter sets for ensemble
params1 = LVector(
    Enz1_Vmax=1.0, Enz1_Keq=10.0, Enz1_K_A_media=1.0, Enz1_K_A=1.0,
    Enz2_Vmax=1.0, Enz2_Keq=1.0, Enz2_K_A=2.0, Enz2_K_B=1.0,
    Enz3_Vmax=1.0, Enz3_Keq=10.0, Enz3_K_B=1.0, Enz3_K_C=1.0,
    Enz4_Vmax=1.0, Enz4_Keq=1.0, Enz4_K_C=1.0, Enz4_K_D=1.0
)

params2 = LVector(
    Enz1_Vmax=2.0, Enz1_Keq=8.0, Enz1_K_A_media=0.5, Enz1_K_A=0.5,
    Enz2_Vmax=1.5, Enz2_Keq=2.0, Enz2_K_A=1.0, Enz2_K_B=2.0,
    Enz3_Vmax=0.8, Enz3_Keq=12.0, Enz3_K_B=1.5, Enz3_K_C=0.8,
    Enz4_Vmax=1.2, Enz4_Keq=0.5, Enz4_K_C=0.8, Enz4_K_D=1.2
)

# Create ensemble problem with multiple parameter sets
ensemble_prob = make_EnsembleProblem(pathway, init_cond, [params1, params2])
ensemble_sol = solve(ensemble_prob, Tsit5(), trajectories=2)

# Plot results
using Plots
plot(ensemble_sol, idxs=(:D), lw=2, alpha=0.7, title="Product D formation", 
     labels=["Set 1" "Set 2"], xlabel="Time", ylabel="[D]")
```

## Uncertainty quantification using Distributions of parameter values

Easily explore parameter space by sampling from statistical distributions:

```julia
using CellMetabolismBase
using DifferentialEquations
using LabelledArrays
using Distributions

# Define a pathway (simplified for clarity)
pathway = MetabolicPathway(
    (:A_media,),
    (
        (:Enz1, (:A_media,), (:A,)),
        (:Enz2, (:A,), (:B,)),
    ),
)

# Define enzyme rate laws
function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz1,(:A_media,),(:A,)}, metabs, params)
    return params.Enz1_Vmax * (metabs.A_media - metabs.A / params.Enz1_Keq)
end

function CellMetabolismBase.enzyme_rate(::Enzyme{:Enz2,(:A,),(:B,)}, metabs, params)
    return params.Enz2_Vmax * metabs.A / (params.Enz2_K_A + metabs.A)
end

# Define initial conditions
init_cond = LVector(A_media=2.0, A=0.0, B=0.0)

# Define parameter distributions
params_dist = LVector(
    Enz1_Vmax = Uniform(0.8, 1.2),    # Uniform distribution
    Enz1_Keq = Normal(10.0, 1.0),     # Normal distribution
    Enz2_Vmax = LogNormal(0.0, 0.2),  # LogNormal distribution
    Enz2_K_A = Uniform(0.5, 1.5)
)

# Create ensemble problem with parameter sampling
ensemble_prob = make_EnsembleProblem(pathway, init_cond, params_dist; n_bootstraps=100)
ensemble_sol = solve(ensemble_prob, Tsit5(), trajectories=100)

# Analyze results
times = 0:0.1:10
b_values = [ensemble_sol[i](times)[3,:] for i in 1:100]  # Get all B values at each time
b_mean = mean(b_values)
b_std = std(b_values)

# Plot results with confidence interval
using Plots
plot(times, b_mean, ribbon=b_std, fillalpha=0.3, 
     title="Mean B production with uncertainty", 
     label="Mean", legend=:bottomright, xlabel="Time", ylabel="[B]")
```

See the [documentation](https://DenisTitovLab.github.io/CellMetabolismBase.jl/stable/) for more detailed examples and API reference.