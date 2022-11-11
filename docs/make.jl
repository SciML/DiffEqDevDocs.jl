using DiffEqDevTools
using OrdinaryDiffEq
using Documenter

DocMeta.setdocmeta!(DiffEqDevTools, :DocTestSetup, :(using DiffEqDevTools); recursive = true)

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "SciML Developer Documentation",
         authors = "Chris Rackauckas",
         modules = [DiffEqDevTools,OrdinaryDiffEq],
         clean = true, doctest = false,
         strict = [
             :doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # Other available options are
             # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
         ],
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://docs.sciml.ai/DiffEqDevDocs/stable/"),
         pages = pages)

deploydocs(;
           repo = "github.com/SciML/DiffEqDevDocs.jl",
           devbranch = "master")

