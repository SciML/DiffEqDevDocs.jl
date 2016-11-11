using Documenter,DiffEqDevTools,DiffEqBase,FiniteElementDiffEq,
      DiffEqProblemLibrary, StokesDiffEq

makedocs(modules=[DiffEqDevTools,DiffEqBase,FiniteElementDiffEq,
                  StokesDiffEq,OrdinaryDiffEq,DiffEqProblemLibrary],
         doctest=false, clean=true,
         format =:html,
         sitename="DifferentialEquations.jl",
         authors="Chris Rackauckas",
         pages = Any[
         "Home" => "index.md",
         "Contributor Guide" => Any[
           "contributing/contributors_guide.md",
           "contributing/adding_algorithms.md",
           "contributing/defining_problems.md",
           "contributing/diffeq_internals.md"
         ],
         "Internal Documentation" => Any[
           "internals/fem_tools.md",
           "internals/extras.md",
           "internals/solver_helpers.md",
           "internals/notes_on_algorithms.md",
           "internals/tableaus.md"
         ]
         ])

deploydocs(
   repo = "github.com/JuliaDiffEq/DiffEqDevDocs.jl.git",
   target = "build",
   osname = "linux",
   julia = "0.5",
   deps = nothing,
   make = nothing)
