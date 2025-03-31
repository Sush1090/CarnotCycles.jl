module CarnotCycles

# Write your package code here.

using ModelingToolkit, CoolProp, DifferentialEquations, Plots, Roots
using Clapeyron
import Clapeyron.PH
import Clapeyron.PS
import ModelingToolkit:scalarize

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

# """
# Storages
# """

include("Components/ThermalStorages/packedbed.jl")


# """
# Plotting
# """

include("Plotting/PhasePlot.jl")
end
