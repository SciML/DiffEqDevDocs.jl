using Documenter,DiffEqDevTools,FiniteElementDiffEq,OrdinaryDiffEq

makedocs(modules=[DiffEqDevTools,FiniteElementDiffEq,OrdinaryDiffEq],
         doctest=false, clean=true,
         format =:html,
         sitename="DiffEq Developer Documentation",
         authors="Chris Rackauckas",
         pages = Any[
         "Home" => "index.md",
         "Contributor Guide" => Any[
           "contributing/ecosystem_overview.md",
           "contributing/adding_packages.md",
           "contributing/adding_algorithms.md",
           "contributing/defining_problems.md",
           "contributing/diffeq_internals.md",
           "contributing/parameters.md",
           "contributing/type_traits.md"
         ],
         "Algorithm Development Tools" => Any[
           "alg_dev/test_problems.md",
           "alg_dev/convergence.md",
           "alg_dev/benchmarks.md"
         ],
         "Internal Documentation" => Any[
           "internals/fem_tools.md",
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
