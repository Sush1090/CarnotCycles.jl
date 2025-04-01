push!(LOAD_PATH,"../src/")
using Documenter,CarnotCycles
makedocs(sitename = "CarnotCycles.jl",
    modules  = [CarnotCycles],
    format = Documenter.HTML(
    # Use clean URLs, unless built as a "local" build
    canonical = "https://sush1090.github.io/CarnotCycles.jl/dev",
),
    pages=[
                "Home" => "index.md"
               ])


deploydocs(;
               repo="github.com/Sush1090/CarnotCycles.jl.git",
           )
           