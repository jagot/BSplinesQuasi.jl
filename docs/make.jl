using Documenter
using BSplinesQuasi

include("plots.jl")

makedocs(
    modules = [BSplinesQuasi],
    sitename = "BSplinesQuasi",
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
    ],
    format = Documenter.HTML(assets = ["assets/latex.js"]),
)

deploydocs(repo = "github.com/jagot/BSplinesQuasi.jl.git")
