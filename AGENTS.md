# Repository Guidelines

## Project Structure & Module Organization
CellMetabolismBase is a Julia package. Core source lives in `src/`, with the module entry point `src/CellMetabolismBase.jl` re-exporting utilities such as `make_ODEProblem` and `make_EnsembleProblem`. Tests are in `test/`, organized as `@testitem` suites per feature (e.g. `tests_make_ODEProblem.jl`, `tests_validate_MetabolicPathway.jl`). Documentation sources sit in `docs/src/`, built via Documenter, and image assets used in the README live under `assets/`. The top-level `Project.toml` defines the package environment; keep `Manifest.toml` in sync when updating dependencies.

## Build, Test, and Development Commands
Run `julia --project=. -e 'using Pkg; Pkg.instantiate()'` after cloning to install dependencies. Execute `julia --project=. -e 'using Pkg; Pkg.test()'` to run the full TestItemRunner-driven suite, including Aqua, JET, and CellMetabolism compatibility checks. For quick iteration you can target a single test item with `julia --project=. -e 'using TestItemRunner; runitems("make_ODEProblem")'`. Rebuild the documentation with `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`.

## Coding Style & Naming Conventions
Follow the Julia style guide: four-space indentation, a 92-character soft limit, and trailing commas on multiline literals. Exported types (e.g. `MetabolicPathway`) use CamelCase; functions and variables use snake_case or `make_`-prefixed verbs (`make_ODEProblem`). Provide docstrings for public APIs and prefer explicit module imports in tests (`using LabelledArrays, BenchmarkTools`). Format using `JuliaFormatter.format(".", indent=4)` when touching multiple files, and run Aqua locally if you change API surface.

## Testing Guidelines
Group related assertions inside `@testitem` blocks, sharing common fixtures via the `setup=[TestMetabolicPathway]` helper defined in `test/define_TestMetabolicPathway.jl`. Aim to cover new enzyme manipulation logic with numerical as well as allocation checks (see `tests_make_ODEProblem.jl` for patterns). Tests may depend on extras (`BenchmarkTools`, `OrdinaryDiffEqFIRK`), so ensure they import lazily within each test item. CI publishes coverage to Codecov; keep branch coverage reasonable and prefer deterministic random seeds.

## Commit & Pull Request Guidelines
Write commits in the imperative mood with concise subjects (`implement disequilibrium_ratios()`) and bundle related changes per commit. Reference issues or PRs with `(#NN)` when applicable. Before opening a pull request, run the full test suite and doc build, summarize functional changes, call out any API breaks, and include reproduction snippets or screenshots for user-visible updates.
