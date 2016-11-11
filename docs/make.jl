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
         "Internal Documentation" => Any[
           "internals/contributors_guide.md",
           "internals/fem_tools.md",
           "internals/extras.md",
           "internals/notes_on_algorithms.md",
           "internals/function_index.md"
         ]
         ])

deploydocs(
   repo = "github.com/JuliaDiffEq/DiffEqDevDocs.jl.git",
   target = "build",
   osname = "linux",
   julia = "0.5",
   deps = nothing,
   make = nothing)
