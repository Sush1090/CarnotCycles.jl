@component function Pump(;name,fluid = set_fluid)
    if fluid isa AbstractString
        return PumpCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return PumpClapeyron(name=name,fluid=set_fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
end

@component function PumpCoolProp(;name,fluid=set_fluid)
    
end


@component function PumpClapeyron(;name,fluid=set_fluid)
    
end