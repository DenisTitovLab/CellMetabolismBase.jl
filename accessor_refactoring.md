# Accessor Refactoring Plan

## Objectives
- Reuse `_generate_Enzymes` to drive both single-enzyme and pathway-wide accessor APIs.
- Provide consistent single-enzyme helpers (`enzyme_name`, `substrates_name`, etc.) that back the pathway aggregators.
- Align `disequilibrium_ratio`/`disequilibrium_ratios` with the same pattern and keep parameter validation behaviour.
- Cover the new accessors with dedicated tests (both singular and aggregate forms) without regressing current benchmarks.

## Proposed Steps
1. **Survey & Prep**
   - Catalogue the existing accessor implementations in `src/MetabolicPathway_and_Enzyme_types_and_related_functions.jl` and `src/make_ODEProblem.jl`.
   - Note which exports and test expectations need to change when the new helpers are introduced.
2. **Implement Single-Enzyme Accessors**
   - Define `enzyme_name`, `substrates_name`, `products_name`, `activators_name`, `inhibitors_name`, and `disequilibrium_ratio` for `Enzyme` objects.
   - Keep argument validation identical to today (especially for equilibrium constants) and add lightweight comments where logic is non-obvious.
3. **Refactor Pathway-Wide Functions**
   - Rewrite the aggregate accessors in `MetabolicPathway_and_Enzyme_types_and_related_functions.jl` to call `_generate_Enzymes` and map over the single-enzyme helpers.
   - Remove and/or adapt the current generated implementations so behaviour stays consistent (ordered tuples, zero allocations).
4. **Update Exports & Tests**
   - Export any newly public helper functions from `src/CellMetabolismBase.jl`.
   - Extend `test/tests_MetabolicPathway_and_Enzyme_types_and_related_functions.jl` to exercise the new single-enzyme APIs and ensure aggregates delegate correctly.
   - Adjust or add benchmarking checks if needed to confirm the refactor keeps performance characteristics.
5. **Verify**
   - Run the relevant test targets (at minimum the accessor testset) and lint for typos or docstring inaccuracies.
   - Stage changes after review and prepare for follow-up feedback before finalizing.

## Progress Log
- [x] 2025-10-15 — Confirmed existing accessor functions and located `_generate_Enzymes` usage in `src/make_ODEProblem.jl` and associated tests.
- [x] 2025-10-15 — Added single-enzyme accessor helpers, centralized `_generate_Enzymes`, and refactored pathway accessors plus disequilibrium ratios to delegate through them; updated exports ready for testing.
- [x] 2025-10-15 — Detected `/usr/local/bin/julia` (v1.10) taking precedence over juliaup shim; reran test suite with `~/.juliaup/bin/julia` (v1.12) and confirmed all tests pass.
- [x] 2025-10-15 — Renamed `reactants_names` accessor to `reactant_names` across source, exports, and tests.
- [x] 2025-10-15 — Renamed pathway aggregate accessors to `substrates_names`, `products_names`, `activators_names`, and `inhibitors_names`; updated exports, validations, and tests to match.
