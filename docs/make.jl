using Documenter, CarnotCycles

makedocs(
sitename="CarnotCycles.jl",
modules = [CarnotCycles],
doctest = true,
pages = [
    "Home" => "index.md",
    "Tutorial" => "Tutorial.md",
    "Cycles Modeling" => "Cycles.md",
    "Cycle Optimization" => "Optimization.md",
    "References" => "reference.md"
]
)

deploydocs(repo = "https://github.com/Sush1090/CarnotCycles.jl",push_preview = true)