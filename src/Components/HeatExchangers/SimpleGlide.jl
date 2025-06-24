"""
This checks if the temperature profile inside the evaporator violates physics or not.

    i.e. It will see for a evaporator that if the temperature of the working fluid was always less than the temperature of the heat transfer fluid.
    if it is feasible then it will return `true` else `false`
"""
function is_feasible_evaporator(T_wf::AbstractVector{T},T_htf::AbstractVector{T}) where T
    n = length(T_wf)
    @assert length(T_htf) == n
    bit_vec = T_htf .> T_wf
    val = prod(bit_vec)
end
@register_symbolic is_feasible_evaporator(T_wf::AbstractVector,T_htf::AbstractVector)

"""
This checks if the temperature profile inside the condensor violates physics or not.

    i.e. It will see for a condensor that if the temperature of the working fluid was always more than the temperature of the heat transfer fluid.
    if it is feasible then it will return `true` else `false`
"""
function is_feasible_condensor(T_wf::AbstractVector,T_htf::AbstractVector)
    n = length(T_wf)
    @assert length(T_htf) == n
    bit_vec = T_htf .< T_wf
    val = prod(bit_vec)
end
@register_symbolic is_feasible_condensor(T_wf::AbstractVector,T_htf::AbstractVector)


"""
A simple evaporator where the HTF inlet and outlet temperature is passed as a parameter.
Has a variable `is_feas` which checks if the fluid passed through violates temperature profile condition or not.
"""
@component function SimpleEvaporatorGlide(;name,fluid=set_fluid,N = 10)
    if fluid isa AbstractString
        return SimpleEvaporatorGlideCoolProp(name=name,fluid=fluid,N = N)
    end
    if fluid isa EoSModel
        return SimpleEvaporatorGlideClapeyron(name=name,fluid=fluid,N = N)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end 
end
export SimpleEvaporatorGlide

@component function SimpleEvaporatorGlideClapeyron(;name,fluid,N)
    @named heatport = HeatPort()
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)

    vars = @variables begin
        is_feas(t)
        T_internal(t)[1:N]
        T_internal_htf(t)[1:N]
        h_internal(t)[1:N]


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
        ΔT_sh(t), [description = "Superheat temperature (K)",bounds = (1e-3,Inf)]
        T_htf_in(t), [description = "Inlet temperature of HTF"]
        T_htf_out(t), [description = "Outlet temperature of HTF"]
    end
    eqs = [
        h_in ~ inport.h
        p_in ~ inport.p
        mdot_in ~ inport.mdot
        x_in ~ inport.x
        z_in ~ mass_to_moles(fluid,x_in,mdot_in)
        s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
        T_in ~ ph_temperature(fluid,p_in,h_in,z_in)
        ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)

        p_crit ~ CriticalPressure(fluid,z_in)
        T_crit ~ CriticalTemperature(fluid,z_in)
        T_bubble ~ Bubble_temperature(fluid,p_in,z_in)
        T_dew ~ Dew_temperature(fluid,p_in,z_in)
        T_sat ~ Base.ifelse(inport.p>=p_crit,T_crit,T_dew)

        T_out ~ ΔT_sh+T_sat
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
        heatport.T_in ~ T_htf_in
        heatport.T_out ~ T_htf_out
        heatport.Q ~ Qdot

        [T_internal_htf[i] ~ T_htf_in - (i-1)/(N-1)*(T_htf_in - T_htf_out) for i = 1:N]
        [h_internal[i] ~ h_out - (i-1)/(N-1)*(h_out -h_in) for i = 1:N]
        [T_internal[i] ~ ph_temperature_qp(fluid,p_in,h_internal[i],z_in) for i = 1:N]
        is_feas ~ is_feasible_evaporator(T_internal,T_internal_htf)
    ]
    compose(System(eqs, t, vars, para;name), inport, outport,heatport)
end


@component function SimpleEvaporatorGlideCoolProp(;name,fluid,N)
    @named heatport = HeatPort()
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)

    vars = @variables begin
        is_feas(t)
        T_internal(t)[1:N]
        T_internal_htf(t)[1:N]
        h_internal(t)[1:N]

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
        T_htf_in(t), [description = "Inlet temperature of HTF"]
        T_htf_out(t), [description = "Outlet temperature of HTF"]
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

        T_out ~ T_sat+ΔT_sh
        p_out ~ p_in
        h_out ~ PropsSI("H","T",T_out,"P",p_out,fluid)
        s_out ~ PropsSI("S","H",h_out,"P",p_out,fluid)
        ρ_out ~ PropsSI("D","H",h_out,"P",p_out,fluid)
        mdot_out ~ mdot_in

        outport.h ~ h_out
        outport.p ~ p_out
        outport.mdot ~ mdot_out

        Qdot ~ mdot_out*(h_out - h_in)
        heatport.T_in ~ T_htf_in
        heatport.T_out ~ T_htf_out
        heatport.Q ~ Qdot

        [T_internal_htf[i] ~ T_htf_in - (i-1)/(N-1)*(T_htf_in - T_htf_out) for i = 1:N]
        [h_internal[i] ~ h_out - (i-1)/(N-1)*(h_out -h_in) for i = 1:N]
        [T_internal[i] ~ PropsSI("T","H",h_internal[i],"P",p_in,fluid) for i = 1:N]
        is_feas ~ is_feasible_evaporator(T_internal,T_internal_htf)
    ]
    compose(System(eqs, t, vars, para;name), inport, outport,heatport)
end


"""
A simple condensor where the HTF inlet and outlet temperature is passed as a parameter.
Has a variable `is_feas` which checks if the fluid passed through violates temperature profile condition or not.
"""
@component function SimpleCondensorGlide(;name,fluid=set_fluid,N = 10)
    if fluid isa AbstractString
        return SimpleCondensorGlideCoolProp(name=name,fluid=fluid,N = N)
    end
    if fluid isa EoSModel
        return SimpleCondensorGlideClapeyron(name=name,fluid=fluid,N = N)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end 
end


@component function SimpleCondensorGlideClapeyron(;name,fluid,N)
    @named heatport = HeatPort()
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)

    vars = @variables begin
        is_feas(t)
        T_internal(t)[1:N]
        T_internal_htf(t)[1:N]
        h_internal(t)[1:N]


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
        T_htf_in(t), [description = "Inlet temperature of HTF"]
        T_htf_out(t), [description = "Outlet temperature of HTF"]
    end
    eqs = [
        h_in ~ inport.h
        p_in ~ inport.p
        mdot_in ~ inport.mdot
        x_in ~ inport.x
        z_in ~ mass_to_moles(fluid,x_in,mdot_in)
        s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
        T_in ~ ph_temperature(fluid,p_in,h_in,z_in)
        ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)

        p_crit ~ CriticalPressure(fluid,z_in)
        T_crit ~ CriticalTemperature(fluid,z_in)
        T_bubble ~ Bubble_temperature(fluid,p_in,z_in)
        T_dew ~ Dew_temperature(fluid,p_in,z_in)
        T_sat ~ Base.ifelse(inport.p>=p_crit,T_crit,T_bubble)

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
        heatport.T_in ~ T_htf_in
        heatport.T_out ~ T_htf_out
        heatport.Q ~ Qdot

        [T_internal_htf[i] ~ T_htf_in - (i-1)/(N-1)*(T_htf_in - T_htf_out) for i = 1:N]
        [h_internal[i] ~ h_out - (i-1)/(N-1)*(h_out -h_in) for i = 1:N]
        [T_internal[i] ~ ph_temperature_qp(fluid,p_in,h_internal[i],z_in) for i = 1:N]
        is_feas ~ is_feasible_condensor(T_internal,T_internal_htf)
    ]
    compose(System(eqs, t, vars, para;name), inport, outport,heatport)
end


@component function SimpleCondensorGlideCoolProp(;name,fluid,N)
    @named heatport = HeatPort()
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)

    vars = @variables begin
        is_feas(t)
        T_internal(t)[1:N]
        T_internal_htf(t)[1:N]
        h_internal(t)[1:N]

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
        T_htf_in(t), [description = "Inlet temperature of HTF"]
        T_htf_out(t), [description = "Outlet temperature of HTF"]
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

        T_out ~ T_sat - ΔT_sc
        p_out ~ p_in
        h_out ~ PropsSI("H","T",T_out,"P",p_out,fluid)
        s_out ~ PropsSI("S","H",h_out,"P",p_out,fluid)
        ρ_out ~ PropsSI("D","H",h_out,"P",p_out,fluid)
        mdot_out ~ mdot_in

        outport.h ~ h_out
        outport.p ~ p_out
        outport.mdot ~ mdot_out

        Qdot ~ mdot_out*(h_out - h_in)
        heatport.T_in ~ T_htf_in
        heatport.T_out ~ T_htf_out
        heatport.Q ~ Qdot

        [T_internal_htf[i] ~ T_htf_in - (i-1)/(N-1)*(T_htf_in - T_htf_out) for i = 1:N]
        [h_internal[i] ~ h_out - (i-1)/(N-1)*(h_out -h_in) for i = 1:N]
        [T_internal[i] ~ PropsSI("T","H",h_internal[i],"P",p_in,fluid) for i = 1:N]
        is_feas ~ is_feasible_condensor(T_internal,T_internal_htf)
    ]
    compose(System(eqs, t, vars, para;name), inport, outport,heatport)
end

export SimpleCondensorGlide