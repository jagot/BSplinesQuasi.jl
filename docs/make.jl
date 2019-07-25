using Documenter
using BSplinesQuasi

include("plots.jl")

makedocs(
    modules = [BSplinesQuasi],
    sitename = "BSplinesQuasi",
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Usage" => [
            "Basis creation" => "usage.md",
            "Splines" => [
                "Spline creation & evaluation" => "splines.md",
                "Function approximation" => "function_approximation.md",
            ],
            "Examples" => [
                "Differentiating functions" => "differentiation.md",
            ]
        ],
    ],
    format = Documenter.HTML(assets = ["assets/latex.js"]),
)

deploydocs(repo = "github.com/jagot/BSplinesQuasi.jl.git")
