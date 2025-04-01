push!(LOAD_PATH,"../src/")
using Documenter,CarnotCycles
makedocs(sitename = "CarnotCycles.jl",
    modules  = [CarnotCycles],
    pages=[
                "Home" => "index.md"
               ])


deploydocs(;
               repo="github.com/Sush1090/CarnotCycles.jl.git",
           )
           