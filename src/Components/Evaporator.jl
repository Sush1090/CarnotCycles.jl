

@component function SimpleEvaporator(;name,fluid=set_fluid)
    if fluid isa AbstractString
        return SimpleEvaporatorCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return SimpleEvaporatorClapeyron(name=name,fluid=fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end 
end


function SimpleEvaporatorCoolProp(;name,fluid)
    @named heatport = HeatPort()
    @named inport = CoolantPort()
    @named outport = CoolantPort()

    vars = @variables begin
        s_in(t)
        h_in(t)
        T_in(t)
        p_in(t)
        ρ_in(t)
        mdot_in(t)

        s_out(t)
        h_out(t)
        T_out(t)
        p_out(t)
        ρ_out(t)
        mdot_out(t)

        p_crit(t)
        T_sat(t)
        T_out(t)

        Qdot(t)
    end
    para = @parameters begin
        ΔT_sh(t), [description = "Superheat temperature (K)",bounds = (1e-3,Inf)]
    end
    eqs = [
        h_in ~ inport.h
        p_in ~ inport.p
        mdot_in ~ inport.mdot
        s_in ~ PropsSI("S","H",h_in,"P",p_in,fluid)
        T_in ~ PropsSI("T","H",h_in,"P",p_in,fluid)
        ρ_in ~ PropsSI("D","H",h_in,"P",p_in,fluid)

        p_crit ~ CritPropSI("PCRIT",fluid)
        T_sat ~ Base.ifelse(inport.p>=p_crit,CritPropSI("TCRIT",fluid),PropsSI("T","Q",1,"P",inport.p,fluid))

        T_out ~ ΔT_sh+T_sat
        p_out ~ p_in
        h_out ~ PropsSI("H","T",T_out,"P",p_out,fluid)
        s_out ~ PropsSI("S","H",h_out,"P",p_out,fluid)
        ρ_out ~ PropsSI("D","H",h_out,"P",p_out,fluid)
        mdot_out ~ mdot_in

        outport.h ~ h_out
        outport.p ~ p_out
        outport.mdot ~ mdot_out

        Qdot ~ mdot_out*(h_out - h_in)
        heatport.T_in ~ T_in
        heatport.T_out ~ T_out
        heatpor.Q ~ Qdot
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport,heatport)
    
end


