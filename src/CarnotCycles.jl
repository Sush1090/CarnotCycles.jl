module CarnotCycles

# Write your package code here.

using ModelingToolkit, CoolProp, DifferentialEquations, Plots
using Clapeyron
import Clapeyron.PH
import Clapeyron.PS
import ModelingToolkit:scalarize

include("ClapeyronRegistrations.jl")

include("Utils.jl")
include("Processes.jl")
include("ThemodynamicClaculations.jl")

include("Components/compressor.jl")
include("Components/expander.jl")
include("Components/HeatExchangers.jl")
include("Components/valve.jl")
include("Components/recuperator.jl")

include("Plotting/PhasePlot.jl")
end
