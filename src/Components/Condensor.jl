

@component function SimpleCondensor(;name,fluid=set_fluid)
    if fluid isa AbstractString
        return SimpleCondensorCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return SimpleCondensorClapeyron(name=name,fluid=fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end 
end


function SimpleCondensorCoolProp(;name,fluid)
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
        ΔT_sc(t), [description = "Subcool temperature (K)",bounds = (1e-3,Inf)]
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

        T_out ~ T_sat-ΔT_sc
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
        heatport.Q ~ Qdot
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport,heatport)
    
end


function SimpleCondensorClapeyron(;name,fluid)
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
        x_in(t)
        z_in(t)

        s_out(t)
        h_out(t)
        T_out(t)
        p_out(t)
        ρ_out(t)
        mdot_out(t)
        x_out(t)
        z_out(t)

        p_crit(t)
        T_sat(t)
        T_out(t)
        T_crit(t)
        T_bubble(t)
        T_dew(t)

        Qdot(t)
    end
    para = @parameters begin
        ΔT_sc(t), [description = "Subcool temperature (K)",bounds = (1e-3,Inf)]
    end
    eqs = [
        h_in ~ inport.h
        p_in ~ inport.p
        mdot_in ~ inport.mdot
        x_in ~ inport.x
        z_in ~ mass_to_moles(fluid,[x_in,1-x_in],mdot_in)
        s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
        T_in ~ ph_temperature(fluid,p_in,h_in,z_in)
        ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)

        p_crit ~ CriticalPressure(fluid,z_in)
        T_crit ~ CriticalTemperature(fluid,z_in)
        T_bubble ~ Bubble_temperature(fluid,p_in,z_in)
        T_dew ~ Dew_temperature(fluid,p_in,z_in)
        T_sat ~ Base.ifelse(inport.p>=p_crit,T_crit,T_dew)

        T_out ~ T_sat - ΔT_sc
        p_out ~ p_in
        z_out ~ z_in
        x_out ~ x_in
        h_out ~ pt_enthalpy(fluid,p_out,T_out,z_out)
        s_out ~ pt_entropy(fluid,p_out,T_out,z_out)
        ρ_out ~ pt_mass_density(fluid,p_out,T_out,z_out)
        mdot_out ~ mdot_in

        outport.h ~ h_out
        outport.p ~ p_out
        outport.mdot ~ mdot_out
        outport.x ~ x_out

        Qdot ~ h_out - h_in
        heatport.T_in ~ T_in
        heatport.T_out ~ T_out
        heatport.Q ~ Qdot
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport,heatport)
    
end

export SimpleCondensor