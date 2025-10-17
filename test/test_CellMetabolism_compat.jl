# @testitem "CellMetabolism glycolysis validates" begin
#     using CellMetabolismBase
#     using CellMetabolism

#     pathway = CellMetabolism.glycolysis_pathway
#     init = deepcopy(CellMetabolism.glycolysis_init_conc)
#     params = deepcopy(CellMetabolism.glycolysis_params)

#     @test isnothing(CellMetabolismBase.validate_MetabolicPathway(pathway, init, params))
# end

# @testitem "CellMetabolism glycolysis make_ODEProblem" begin
#     using CellMetabolism
#     using OrdinaryDiffEqFIRK
#     using SciMLBase

#     pathway = CellMetabolism.glycolysis_pathway
#     init = deepcopy(CellMetabolism.glycolysis_init_conc)
#     params = deepcopy(CellMetabolism.glycolysis_params)
#     tspan = (0.0, 60.0)

#     prob = CellMetabolismBase.make_ODEProblem(pathway, init, tspan, params)
#     @test prob.tspan == tspan
#     @test prob.u0 == init
#     @test prob.p == params

#     sol = solve(prob, RadauIIA9(); abstol = 1e-15, reltol = 1e-8, save_everystep = false)
#     @test sol.retcode === SciMLBase.ReturnCode.Success
# end

# @testitem "CellMetabolism glycolysis make_EnsembleProblem" begin
#     using CellMetabolism
#     using OrdinaryDiffEqFIRK
#     using SciMLBase

#     pathway = CellMetabolism.glycolysis_pathway
#     base_init = deepcopy(CellMetabolism.glycolysis_init_conc)
#     base_params = deepcopy(CellMetabolism.glycolysis_params)
#     tspan = (0.0, 30.0)

#     init_ensemble = [scale .* base_init for scale in (0.95, 1.05)]
#     params_ensemble = [scale .* base_params for scale in (0.9, 1.1)]

#     ensemble_prob = CellMetabolismBase.make_EnsembleProblem(pathway, init_ensemble, params_ensemble; tspan = tspan)
#     sols = solve(
#         ensemble_prob,
#         RadauIIA9();
#         trajectories = length(init_ensemble),
#         save_everystep = false,
#         abstol = 1e-15,
#         reltol = 1e-8,
#     )

#     @test length(sols) == length(init_ensemble)
#     @test all(sol.retcode == SciMLBase.ReturnCode.Success for sol in sols)
# end
