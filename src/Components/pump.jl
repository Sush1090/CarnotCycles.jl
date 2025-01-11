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
    @named inport =CoolantPort()
    @named outport = CoolantPort()
    vars = @variables begin
        LiquidPhase(t)
        s_in(t)
        T_in(t)
        p_in(t)
        h_in(t)
        ρ_in(t)
        mdot_in(t)
        
        s_out(t)
        T_out(t)
        p_out(t)
        h_out(t)
        ρ_out(t)
        mdot_out(t)
        P(t)
    end
    para = @parameters begin
        η(t),[description = "Isentropic Efficiency"]
        πc(t), [description = "Pressure Ratio"]
    end
    eqs = [
        h_in ~ inport.h
        p_in ~ inport.p
        mdot_in ~ inport.mdot
        s_in ~ PropsSI("S","H",h_in,"P",p_in,fluid)
        T_in ~ PropsSI("T","H",h_in,"P",p_in,fluid)
        ρ_in ~ PropsSI("D","H",h_in,"P",p_in,fluid)
        LiquidPhase ~ CoolPropLiquidPhaseCheck(fluid,h,p)

        mdot_out ~ mdot_in
        p_out ~ p_in*πc
        h_out ~ IsentropicCompression(πc, h_in, p_in,fluid,η)
        s_out ~ PropsSI("S","H",h_out,"P",p_out,fluid)
        T_out ~ PropsSI("T","H",h_out,"P",p_out,fluid)
        ρ_out ~ PropsSI("D","H",h_out,"P",p_out,fluid)
        
        P ~ mdot_in*(h_out - h_in)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end


@component function PumpClapeyron(;name,fluid=set_fluid)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    vars = @variables begin
        LiquidPhase(t)
        s_in(t)
        T_in(t)
        p_in(t)
        h_in(t)
        ρ_in(t)
        z_in(t)
        x_in(t)
        mdot_in(t)

        s_out(t)
        T_out(t)
        p_out(t)
        h_out(t)
        ρ_out(t)
        z_out(t)
        x_out(t)
        mdot_out(t)

        P(t)
    end
    para = @parameters begin
        πc(t), [description = "Pressure ratio (-)"]
        η(t), [description = "Isentropic Efficiency (-)"]
    end

    eqs = [
            LiquidPhase ~ LiquidPhaseChecker(fluid,p_in,h_in,z_in) 
            mdot_in ~ inport.mdot
            x_in ~ inport.x
            z_in ~ mass_to_moles(fluid,x_in,mdot_in)
            p_in ~ inport.p
            h_in ~ inport.h
            s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
            T_in ~ ph_temperature(fluid,p_in,h_in,z_in)
            ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)

            x_out ~ x_in
            z_out ~ z_in
            p_out ~ p_in*πc
            h_out ~ IsentropicCompressionClapeyron(πc, h_in, p_in,z_in,fluid,η)
            s_out ~ ph_entropy(fluid,p_out,h_out,z_out)
            T_out ~ ph_temperature(fluid,p_out,h_out,z_out)
            ρ_out ~ ph_mass_density(fluid,p_out,h_out,z_out)
            mdot_out ~ mdot_in

            outport.mdot ~ mdot_out
            outport.x ~ x_out
            outport.p ~ p_out
            outport.h ~ h_out

            P ~ h_out - h_in
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

export Pump