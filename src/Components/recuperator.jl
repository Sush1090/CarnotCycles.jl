

abstract type RecuperatoreType end

struct RecuperatorAverageTemp <:RecuperatoreType 

end

struct RecuperatorORC <:RecuperatoreType 

end
export RecuperatorORC

function Recuperator_TemperatureOut(T_gas_in,T_liquid_in,T_sat,ΔT_sat_diff)
    if (T_gas_in - ΔT_sat_diff <= T_sat)
        return T_gas_in;
    end
    if (T_gas_in < T_liquid_in)
        return T_gas_in;
    end
    if (T_gas_in - ΔT_sat_diff <= T_liquid_in)
        return T_gas_in;
    end

    T_gas_out =  T_gas_in - ΔT_sat_diff
    return T_gas_out
end
@register_symbolic Recuperator_TemperatureOut(T_gas_in,T_liquid_in,T_sat,ΔT_sat_diff)
export Recuperator_TemperatureOut

function Recuperator(type::RecuperatorORC;name,ΔT_sat_diff,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @assert ΔT_sat_diff > 1e-3 "Keep ΔT_sat_diff > 1e-3, for safety from saturation curve of CoolProp"
    @named inport_gas = CoolantPort(fluid=fluid)
    @named outport_gas = CoolantPort(fluid=fluid)
    @named inport_liquid = CoolantPort(fluid=fluid)
    @named outport_liquid = CoolantPort(fluid=fluid)

    vars = @variables begin
        Q_dot_transfer(t)

        #gas port
        h_in_gas(t)
        h_out_gas(t)
        p_in_gas(t)
        p_out_gas(t)
        T_in_gas(t)
        T_out_gas(t)
        ρ_out_gas(t)
        ρ_in_gas(t)
        s_in_gas(t)
        s_out_gas(t)

        #two phase var
        T_sat(t)
        #liquid port
        h_in_liquid(t)
        h_out_liquid(t)
        p_in_liquid(t)
        p_out_liquid(t)
        T_in_liquid(t)
        T_out_liquid(t)
        ρ_out_liquid(t)
        ρ_in_liquid(t)
        s_in_liquid(t)
        s_out_liquid(t)
    end

    para = @parameters begin
        
    end

    eqs = [
        h_in_gas ~ inport_gas.h
        p_in_gas ~ inport_gas.p
        outport_gas.mdot ~ abs(inport_gas.mdot)

        h_in_liquid ~ inport_liquid.h
        p_in_liquid ~ inport_liquid.p 
        outport_liquid.mdot ~ abs(inport_liquid.mdot)

        T_sat ~ PropsSI("T","Q",1,"P",inport_gas.p,fluid)
        T_in_gas ~ PropsSI("T","H",inport_gas.h,"P",inport_gas.p,fluid)
        T_in_liquid ~ PropsSI("T","H",inport_liquid.h,"P",inport_liquid.p,fluid)
        T_out_gas ~ Recuperator_TemperatureOut(T_in_gas,T_in_liquid,T_sat,ΔT_sat_diff)

        p_out_gas ~ p_in_gas
        p_out_liquid ~ p_in_liquid

        h_out_gas ~ PropsSI("H","T",T_out_gas,"P",p_out_gas,fluid)
        outport_gas.h ~ h_out_gas
        outport_gas.p ~ p_out_gas


        Q_dot_transfer ~ outport_gas.mdot*(h_out_gas - h_in_gas)

        h_out_liquid ~ h_in_liquid - Q_dot_transfer/outport_liquid.mdot
        outport_liquid.h ~ h_out_liquid
        outport_liquid.p ~ p_out_liquid

        s_in_gas ~ PropsSI("S","H",h_in_gas,"P",p_in_gas,fluid)
        ρ_in_gas ~ PropsSI("D","H",h_in_gas,"P",p_in_gas,fluid)
        s_out_gas ~ PropsSI("S","H",h_out_gas,"P",p_out_gas,fluid)
        ρ_out_gas ~ PropsSI("D","H",h_out_gas,"P",p_out_gas,fluid)

        s_in_liquid ~ PropsSI("S","H",h_in_liquid,"P",p_in_liquid,fluid)
        ρ_in_liquid ~ PropsSI("D","H",h_in_liquid,"P",p_in_liquid,fluid)
        s_out_liquid ~ PropsSI("S","H",h_out_liquid,"P",p_out_liquid,fluid)
        ρ_out_liquid ~ PropsSI("D","H",h_out_liquid,"P",p_out_liquid,fluid)
        T_out_liquid ~ PropsSI("T","H",h_out_liquid,"P",p_out_liquid,fluid)
    ]
    compose(System(eqs, t, vars, para;name), inport_gas,inport_liquid,outport_gas,outport_liquid)
end
export Recuperator