using Documenter
using BSplinesQuasi

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

DocMeta.setdocmeta!(BSplinesQuasi, :DocTestSetup, :(using BSplinesQuasi); recursive=true)
makedocs(
    modules = [BSplinesQuasi],
    sitename = "BSplinesQuasi",
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Usage" => [
            "Basis creation" => "usage.md",
            "Knot sets" => "knot_sets.md",
            "Splines" => [
                "Spline creation & evaluation" => "splines.md",
                "Function approximation" => "function_approximation.md",
            ],
            "Examples" => [
                "Differentiating functions" => "differentiation.md",
                "Ordinary differential equations" => "odes.md",
            ]
        ],
    ],
    format = Documenter.HTML(assets = ["assets/latex.js"]),
)

deploydocs(repo = "github.com/jagot/BSplinesQuasi.jl.git")
