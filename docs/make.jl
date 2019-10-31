using Documenter,DiffEqDevTools,OrdinaryDiffEq

makedocs(modules=[DiffEqDevTools,OrdinaryDiffEq],
         doctest=false, clean=true,
         format = Documenter.HTML(),
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
           "contributing/type_traits.md"
         ],
         "Algorithm Development Tools" => Any[
           "alg_dev/test_problems.md",
           "alg_dev/convergence.md",
           "alg_dev/benchmarks.md"
         ],
         "Internal Documentation" => Any[
           "internals/notes_on_algorithms.md",
           "internals/tableaus.md"
         ]
         ])

deploydocs(
   repo = "github.com/JuliaDiffEq/DiffEqDevDocs.jl.git",
   target = "build",
   osname = "linux",
   julia = "1.1",
   deps = nothing,
   make = nothing)
