

ModelingToolkit.@independent_variables t
D = Differential(t)

PropsSI(out::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString) = CoolProp.PropsSI(out, name1, value1, name2, value2, fluid)
@register_symbolic PropsSI(out::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)

PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString) = CoolProp.PhaseSI(name1, value1, name2, value2, fluid)
@register_symbolic PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)

CritPropSI(property::AbstractString,fluid::AbstractString) = CoolProp.PropsSI(property,fluid)
@register_symbolic CritPropSI(property::AbstractString,fluid::AbstractString)

global set_fluid = nothing
global Nc = nothing

"""
`load_fluid(x::AbstractString)` - fixes fluid for simulation through components using CoolProp as backend.
"""
function load_fluid(x::AbstractString)
    CarnotCycles.set_fluid = x
    CarnotCycles.Nc = 1
    return CarnotCycles.set_fluid
end

"""
`load_fluid(x::Clapeyron.EoSModel)` - fixes fluid for simulation through components using Clapeyron as backend
"""
function load_fluid(x::Clapeyron.EoSModel)
    CarnotCycles.set_fluid = x
    CarnotCycles.Nc = size(x.components,1)
    @assert size(x.components,1) <=2 "For now we support only two component mixtures"
    return CarnotCycles.set_fluid
end

function load_fluid(x)
    throw(error("The type of fluid provided is not supported yet. It has to be `Clapeyron.EoSModel` or `AbstractString`"))
end

export load_fluid


function Show_fluid_details(fluid=set_fluid)
    if fluid isa AbstractString

    end

    if fluid isa EoSModel

    end
end
export Show_fluid_details

"""
`mass_to_moles(model::EoSModel,x,mass)` : convert mass of fluid to number of moles based on the composition of 1st fluid by mass `x`
"""
function mass_to_moles(model::EoSModel,x,mass)
    @assert 0<x<=1
    if size(model.components,1) == 1
        mws = Clapeyron.mw(model);
        x = 1

        M = mass*x;
        moles =  M./mws
        return [moles[1]]
    end
    x_ = [x,1-x]
    mws = Clapeyron.mw(model);
    M = mass*x_;

    moles =  M./mws
    return moles
end
@register_symbolic mass_to_moles(model::EoSModel,x,mass)



"""
`MassSource` -  Initilizes cycle start point. Requires initial enthalpy,pressure and Mass flow rate if `CoolProp` fluid is used else uses enthalpy,pressure, Mass flow rate, and mass fraction of first fluid if `Clapeyron` Fluid is used.
"""
@component function MassSource(;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    if fluid isa AbstractString
        return MassSourceCoolProp(name=name,fluid = fluid)
    end

    if fluid isa EoSModel
        return MassSourceClapeyron(name=name,fluid = fluid)
    end
end

@component function MassSourceCoolProp(;name, fluid = set_fluid)
    @named port = CoolantPort()
    para = @parameters begin
        source_pressure(t)
        source_enthalpy(t)
        source_mdot(t)
    end
    vars = @variables begin
        mdot(t)
        s(t)
        p(t)
        T(t)
        h(t)
        ρ(t)
     end

    eqs = [
        port.mdot ~ source_mdot # Outflow is negative
        port.p ~ source_pressure
        port.h ~ source_enthalpy
        mdot ~ port.mdot
        s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
        p ~ port.p
        T ~ PropsSI("T","H",port.h,"P",port.p,fluid)
        h ~ port.h
        ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name),port)
end

@component function MassSourceClapeyron(;name, fluid = set_fluid,Nc = Nc) 
    @named port = CoolantPort()
    para = @parameters begin
        source_pressure(t), [description = "Pressure at source (Pa)"]
        source_enthalpy(t), [description = "Enthalpy at source (J)"]
        source_mdot(t), [description = "Moles at source (-)"]
        source_x(t), [description = "mass % of 1st component", bounds = (0,1)]
    end
    vars = @variables begin
        mdot(t), [description = "Mass flow rate (g/s)"]
        x(t), [description = "mass % of 1st component"]
        z(t) , [description = "Moles (-)"]
        s(t), [description = "Entropy (J/mol.K)"]
        p(t), [description = "Pressure (Pa)"]
        T(t), [description = "Temperature (K)"]
        h(t), [description = "Enthalpy (J)"]
        ρ(t), [description = "Total Density (kg)"]
     end

    eqs = [
        port.mdot ~ source_mdot
        port.p ~ source_pressure
        port.h ~ source_enthalpy
        port.x ~ source_x
        mdot ~ source_mdot
        x ~ source_x
        z ~ mass_to_moles(fluid,x,mdot)
        s ~ ph_entropy(fluid,p,h,z)
        p ~ port.p
        T ~ ph_temperature(fluid,p,h,z)
        h ~ port.h
        ρ ~ ph_mass_density(fluid,p,h,z)
    ]
    compose(ODESystem(eqs, t, vars, para;name),port)
end


"""
`MassSink`` -  Sets the final port input values to the variables. Use as a sink for the simulation.
"""
@component function MassSink(;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    if fluid isa AbstractString
        return MassSinkCoolProp(name=name,fluid = fluid)
    end

    if fluid isa EoSModel
        return MassSinkClapeyron(name=name,fluid = fluid)
    end
end

function MassSinkCoolProp(;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named    port = CoolantPort()
    para = @parameters begin
        
    end
    vars = @variables begin
        mdot(t)
        s(t)
        p(t)
        T(t)
        h(t)
        ρ(t)
     end

   eqs = [
    port.p ~ p
    port.h ~ h
    mdot ~ port.mdot
    s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
    T ~ PropsSI("T","H",port.h,"P",port.p,fluid)
    ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name),port)
end

@component function MassSinkClapeyron(;name,fluid = set_fluid)
        @named    port = CoolantPort()
        para = @parameters begin
            
        end
        vars = @variables begin
            mdot(t)
            x(t)
            z(t)
            s(t)
            p(t)
            T(t)
            h(t)
            ρ(t)
         end
    
       eqs = [
        port.p ~ p
        port.h ~ h
        x ~ port.x
        mdot ~ port.mdot
        z ~ mass_to_moles(fluid,x,mdot)
        s ~ ph_entropy(fluid,p,h,z)
        T ~ ph_temperature(fluid,p,h,z)
        ρ ~ ph_mass_density(fluid,p,h,z)
       ]
       compose(ODESystem(eqs, t, vars, para;name),port)
end


export CoolantPort,CoolComponent,MassSink,AmbientTemperature,AtmosphericPressure,ComputeSpecificLatentHeat
export MassSource






"""
`CoolPropLiquidPhaseCheck` Checks if state point is Liquid. Passes through an assert statement. 
"""
function CoolPropLiquidPhaseCheck(fluid::AbstractString,h,p)
    phase = PhaseSI("H",h,"P",p,fluid)
    @assert phase == "liquid" || phase == "supercritical_liquid"
    return 1;
end
@register_symbolic CoolPropLiquidPhaseCheck(fluid::AbstractString,h,p)

"""
`CoolPropGasPhaseCheck` Checks if state point is Gas. Passes through an assert statement. 
"""
function CoolPropGasPhaseCheck(fluid::AbstractString,h,p)
    phase = PhaseSI("H",h,"P",p,fluid)
    @assert phase == "gas" || phase == "supercritical"
    return 1;
end
@register_symbolic CoolPropGasPhaseCheck(fluid::AbstractString,h,p)



@component function Heat2Storage(;name)
    @named heatport = HeatPort()
    @named storeport = StoragePort()
    para = @parameters begin
       Cp(t) , [description = "Specific heat of HTF (J/Kg/K)"]
       T_in(t), [description  = "Inlet Temperature of HTF (K)"]
       T_out(t), [description  = "Outlet Temperature of HTF (K)"]
    end
    vars = @variables begin
       mdot(t), [description = "masss flow Rate of HTF (kg/s)"]
    end
    eqs = [
        storeport.mdot ~ mdot
        mdot ~ -heatport.Q/(Cp*(T_out - T_in))
        storeport.T ~ T_out
    ]
    return compose(ODESystem(eqs, t, vars, para;name=name),heatport,storeport)
end

export Heat2Storage


"""
`show_all_states(sol::SteadyStateSolution,system::Vector{ODESystem},names::Vector{String})`

* Arguments:
    1. `sol`    : The solution from the  SteadtStateSystem. 
    2. `system` : The vector of components chosen. The first should be `MassSource` and the last has to be `MassSink`.
    3. `names`  : The names of the components in the system. Pass it as a vector of strings.

* Returns: Prints the state points of the system in a readable format in the terminal.
"""
function show_all_states(sol::SteadyStateSolution,system::Vector{ODESystem},names::Vector{String})

    @assert length(system)-2 == length(names) "Length of system and names must be equal, excluding source and sink"
    system = system[2:end-1] # Exclude source and sink
    if CarnotCycles.set_fluid isa EoSModel
        for i in eachindex(names)
            println("--------------------------------------")
            println("Component: ", names[i])
            println("Temperature in : ", sol[system[i].T_in], " K ", " , " , "Temperature out : ", sol[system[i].T_out], " K")
            println("Pressure in : ", sol[system[i].p_in], " Pa ", " , " , "Pressure out : ", sol[system[i].p_out], " Pa")
            println("Enthalpy : ", sol[system[i].h_in], " J ", " , " , "Enthalpy out : ", sol[system[i].h_out], " J")
            println("Entropy : ", sol[system[i].s_in], " J/K ", " , " , "Entropy out : ", sol[system[i].s_out], " J/K")
            println("Mass Flow Rate : ", sol[system[i].inport.mdot], " g/s ", " , " , "Mass Flow Rate out : ", sol[system[i].outport.mdot], " g/s")
          
            println("Mass fraction : ", sol[system[i].inport.x], " , " , "Mass fraction out : ", sol[system[i].outport.x])
           
            println("--------------------------------------")
        end
    end

    if CarnotCycles.set_fluid isa AbstractString
        for i in eachindex(names)
            println("--------------------------------------")
            println("Component: ", names[i])
            println("Temperature in : ", sol[system[i].T_in], " K ", " , " , "Temperature out : ", sol[system[i].T_out], " K")
            println("Pressure in : ", sol[system[i].p_in], " Pa ", " , " , "Pressure out : ", sol[system[i].p_out], " Pa")
            println("Enthalpy : ", sol[system[i].h_in], " J/kg ", " , " , "Enthalpy out : ", sol[system[i].h_out], " J/kg")
            println("Entropy : ", sol[system[i].s_in], " J/K/kg ", " , " , "Entropy out : ", sol[system[i].s_out], " J/K/kg")
            println("Mass Flow Rate : ", sol[system[i].inport.mdot], " kg/s ", " , " , "Mass Flow Rate out : ", sol[system[i].outport.mdot], " kg/s")           
            println("--------------------------------------")
        end
    end
    
end

export show_all_states
