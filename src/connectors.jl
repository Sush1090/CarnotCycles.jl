# """
# Makes single node at ports. This node is Pressure,Enthalpy and Massflowrate
# """
@connector  function CoolantPort(;name,fluid = set_fluid) 
    if fluid isa EoSModel
        return CoolantPortClapeyron(;name = name)
    end
    if fluid isa AbstractString
        return CoolantPortCoolProp(name = name)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
end

@connector function CoolantPortCoolProp(;name)
    vars = @variables begin 
        p(t),  [description = "Pressure (Pa)",input = true]
        h(t), [description = "Enthalpy (J/kg)",input = true]
        mdot(t), [description = "Mass Flow Rate (kg/s)",input = true] 
    end
    return ODESystem(Equation[], t, vars, [];name=name)
end

@connector function CoolantPortClapeyron(;name,Nc=Nc)
    vars = @variables begin 
        p(t),  [description = "Pressure (Pa)",input = true]
        h(t), [description = "Enthalpy (J/kg)",input = true]
        mdot(t), [description = "Mass Flow Rate (g/s)",input = true] 
        x(t), [description = "mass % of 1st component", input = true]
    end
    return ODESystem(Equation[], t, vars, [];name=name)
end

@connector  function RefPortCoolProp(;name) 
    vars = @variables begin 
        p(t),  [description = "Pressure (Pa)",input = true]
        T(t), [description = "Temperature",input = true]
        mdot(t), [description = "Mass Flow Rate (kg/s)",input = true] 
    end
    ODESystem(Equation[], t, vars, [];name=name)
end

@connector  function RefPortClapeyron(;name) 
    vars = @variables begin 
        p(t),  [description = "Pressure (Pa)",input = true]
        T(t), [description = "Temperature",input = true]
        mdot(t), [description = "Mass Flow Rate (g/s)",input = true] 
        x(t), [description = "mass % of 1st component", input = true]
    end
    ODESystem(Equation[], t, vars, [];name=name)
end

@connector function RefPort(;name,fluid = set_fluid)
    if fluid isa EoSModel
        return RefPortClapeyron(;name = name)
    end
    if fluid isa AbstractString
        return RefPortCoolProp(name = name)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
end

@connector function PowerPort(;name)
    vars = @variables begin 
        P(t),  [description = "Power (W)",input = true]
    end
    ODESystem(Equation[], t, vars, [];name=name)
end


@connector function HeatPort(;name)
    vars = @variables begin 
        Q(t),  [description = "Heat rate (W)",input = true]
        T_in(t), [description = "Inlet Temperature of working fluid",input = true]
        T_out(t), [description = "Outlet Temperature of working fluid",input = true]
    end
    ODESystem(Equation[], t, vars, [];name=name)
end

# """
#  Storage port that connect the storage HTF to the thermal storage
# """
@connector function StoragePort(;name)
    vars = @variables begin
        T(t), [input = true,description ="Temperature of Storage HTF"]
        mdot(t), [input = true,description ="Mass Flow Rate of Storage HTF"]
    end
    ODESystem(Equation[],t,vars,[],name=name)
end

export StoragePort

@connector function AmbientNode(;node)
    vars = @variables begin
        T(t), [input = true,description ="Ambient Temperature"]
    end
    ODESystem(Equation[],t,vars,[],name=name)
end

export PowerPort, AmbientNode, HeatPort