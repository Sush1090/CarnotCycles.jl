using ModelingToolkit

ModelingToolkit.@independent_variables t
D = Differential(t)
const AtmosphericPressure = 101305 #Pa
const AmbientTemperature = 300 #K

PropsSI(out::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString) = CoolProp.PropsSI(out, name1, value1, name2, value2, fluid)
@register_symbolic PropsSI(out::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)

PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString) = CoolProp.PhaseSI(name1, value1, name2, value2, fluid)
@register_symbolic PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)


global set_fluid = nothing

"""
`@load_fluid` - fixes fluid for simulation through components
"""
macro load_fluid(x::AbstractString)
    CarnotCycles.set_fluid = x
    return CarnotCycles.set_fluid
end
export @load_fluid

"""
Collection of state point
"""
struct ThermoState
    mdot
    h
    p
    fluid
    function ThermoState(;mdot,h,p,fluid=set_fluid)
        new(mdot,h,p,fluid)
    end
end
export ThermoState
"""
    Initial state point to `MassSource`
"""
function initialize_state(;T_start=nothing,p_start=nothing,mdot=nothing,fluid = set_fluid)
    if isnothing(fluid) == true
        throw(error("Fluid not selected"))
    end
    if isnothing(T_start)
        throw(error("Initial Temperature not chose"))
    end
    if isnothing(p_start)
        throw(error("Initial Pressure not chose"))
    end
    if isnothing(mdot)
        throw(error("mass flow rate not chose"))
    end
    h_start = PropsSI("H","T",T_start,"P",p_start,fluid)
    return ThermoState(mdot = mdot,h=h_start,p = p_start,fluid=fluid)
end
export initialize_state
"""
Makes single node at ports. This node is Pressure,Enthalpy and Massflowrate
"""
@connector  function CoolantPort(;name) 
    vars = @variables begin 
        p(t),  [description = "Pressure (Pa)",input = true]
        h(t), [description = "Enthalpy (J/kg)",input = true]
        mdot(t), [description = "Mass Flow Rate (kg/s)",input = true] # alternates sign of variable after every component.
    end
    ODESystem(Equation[], t, vars, [];name=name)
end


"""
This is double node at ports. Inport and outport. 
"""
function CoolComponent(;name) 
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    vars = @variables begin
        Δp(t) 
        Δmdot(t) 
        Δh(t) 
    end
    eqs = [
        Δp ~ outport.p - inport.p
        Δh ~ outport.h - inport.h
        Δmdot ~ outport.mdot - inport.mdot
    ]
    compose(ODESystem(eqs, t, vars, [];name=name), inport, outport)
end


"""
Mass source -  Use when the cycle needs a start point. Requires initial enthalpy,pressure and Massflowrate
"""
function MassSource(;name,source_pressure = 101305,source_enthalpy=1e6,source_mdot=5,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named port = CoolantPort()
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

function MassSource(state::ThermoState;name)
    @unpack mdot,h,p,fluid = state
    return MassSource(name = name,source_pressure=p,source_enthalpy=h,source_mdot=mdot,fluid=fluid)
end

"""
Mass sink -  Use when the cycle needs a end point. Sets the final port input values to the variables
"""
function MassSink(;name,fluid = set_fluid) 
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
    p ~ port.p
    T ~ PropsSI("T","H",port.h,"P",port.p,fluid)
    h ~ port.h
    ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name),port)
end
"""
`MassSink(state::ThermoState;name,fluid=set_fluid)`

"""
function MassSink(state::ThermoState;name,fluid=set_fluid)
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
    port.p ~ state.p
    port.h ~ h
    mdot ~ state.mdot
    s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
    p ~ state.p
    T ~ PropsSI("T","H",port.h,"P",port.p,fluid)
    h ~ state.h
    ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name),port)
end

"""
ComputeSpecificLatentHeat: Computes the specific latent heat of the give fluid at a particular varliable value. var1 should not be enthalpy or vapour quality
*Arguments:
-'var1'     : Variable String. Uses CoolProp variable strings
-'value1'   : Value of the variale chosen
-'fluid'    : Fluid name string
"""
function ComputeSpecificLatentHeat(var1::AbstractString,value1,fluid::AbstractString)
    @assert var1 != "Q"
    H_L = PropsSI("H",var1,value1,"Q",0,fluid)
    H_V = PropsSI("H",var1,value1,"Q",1,fluid)
    return H_V - H_L
end
@register_symbolic ComputeSpecificLatentHeat(var1::AbstractString,value1,fluid::AbstractString)

export CoolantPort,CoolComponent,MassSink,AmbientTemperature,AtmosphericPressure,ComputeSpecificLatentHeat


function MassSourceLiquid(;name,fluid=set_fluid,Δp_super = 1e3,source_mdot,source_temperature)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    if isnothing(source_mdot)
        throw(error("mass flow rate at source not chosen"))
    end
    if isnothing(source_temperature)
        throw(error("temperature at source not chosen"))
    end
    if isnothing(Δp_super)
        throw(error("pressure above saturation pressure not chosen at liquid source"))
    end
    @named port = CoolantPort()
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
        port.mdot ~ source_mdot # Outflow is negative
        T ~ source_temperature
        port.p ~ PropsSI("P","T",T,"Q",0,fluid) + Δp_super
        port.h ~ PropsSI("H","P",port.p,"T",T,fluid)
        mdot ~ port.mdot
        s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
        p ~ port.p
        h ~ port.h
        ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name),port)
end


function MassSourceGas(;name,fluid=set_fluid,Δp_sub = 1e3,source_mdot,source_temperature)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    if isnothing(source_mdot)
        throw(error("mass flow rate at source not chosen"))
    end
    if isnothing(source_temperature)
        throw(error("temperature at source not chosen"))
    end
    if isnothing(Δp_sub)
        throw(error("pressure below saturation pressure not chosen at gas source"))
    end
    @named port = CoolantPort()
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
        port.mdot ~ source_mdot # Outflow is negative
        T ~ source_temperature
        port.p ~ PropsSI("P","T",T,"Q",0,fluid) - Δp_sub
        port.h ~ PropsSI("H","P",port.p,"T",T,fluid)
        mdot ~ port.mdot
        s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
        p ~ port.p
        h ~ port.h
        ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name),port)
end



"""
`MassSource(type::Symbol; fluid= set_fluid,Δp_sub = nothing,Δp_super = nothing,source_mdot = nothing,source_temperature = nothing)`

Initilize source with temperature and pressure away from saturation curve. Set type to be a symbol between `:liquid` or `:gas`.
*    Arguments: 
    1. `Δp_sub`             : pressure below saturation curve i.e. for gas as source state.
    2. `Δp_super`           : pressure above saturation curve i.e. for liquid as source state.   
    3. `source_mdot`        : mass flow rate of the system
    4. `source_temperature` : source temperature
    5. `type`               : `:liquid` or `:gas`
"""
function MassSource(type::Symbol; fluid= set_fluid,Δp_sub = nothing,Δp_super = nothing,source_mdot = nothing,source_temperature = nothing,name)
    if type == :liquid
       @named src = MassSourceLiquid(fluid=fluid,Δp_super = Δp_super,source_mdot = source_mdot,source_temperature = source_temperature)
       return src
    elseif type == :gas
        @named src = MassSourceGas(fluid=fluid,Δp_sub = Δp_sub,source_mdot = source_mdot,source_temperature = source_temperature)
        return src
    else
        throw(error("Choose the source type to be `:liquid` or `:gas`"))
    end
end

export MassSourceGas, MassSourceLiquid, MassSource