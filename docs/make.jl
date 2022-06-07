using Documenter,DiffEqDevTools,OrdinaryDiffEq

include("pages.jl")

makedocs(modules=[DiffEqDevTools,OrdinaryDiffEq],
         doctest=false, clean=true,
         format = Documenter.HTML(),
         sitename="SciML Developer Documentation",
         authors="Chris Rackauckas",
         pages = pages)

deploydocs(
   repo = "github.com/SciML/DiffEqDevDocs.jl.git",
)
