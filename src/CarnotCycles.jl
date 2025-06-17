module CarnotCycles

# Write your package code here.

using ModelingToolkit, CoolProp, Roots
using Clapeyron
import Clapeyron.PH
import Clapeyron.PS
import ModelingToolkit:scalarize
import Plots: plot
import SteadyStateDiffEq:SteadyStateSolution
# """
# Symbolic Registrations
# """

include("ClapeyronRegistrations.jl")

# """
# Utilities
# """

include("connectors.jl")
include("Utils.jl")
include("Processes.jl")
include("ThemodynamicClaculations.jl")

# """
# Basic Components
# """

include("Components/compressor.jl")
include("Components/expander.jl")
include("Components/HeatSources.jl")
include("Components/valve.jl")
include("Components/recuperator.jl")
include("Components/Evaporator.jl")
include("Components/Condensor.jl")
include("Components/pump.jl")

# """
# Simple HEX models
# """

include("Components/HeatExchangers/SimpleGlide.jl")
include("Components/HeatExchangers/HeatExchangerBase.jl")
include("Components/HeatExchangers/eNTU.jl")


# """
# Storages
# """

include("Components/ThermalStorages/packedbed.jl")


# """
# Pipes
# """
include("Components/pipes.jl")

# """
# Plotting
# """

include("Plotting/PhasePlot.jl")
end
