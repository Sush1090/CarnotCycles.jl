
"""
`PipeCoolProp(fluid::AbstractString = set_fluid;name)`

pressure drop across pipe using Darcy-Weisbach equation
"""
@component function PipeCoolProp(fluid::AbstractString = set_fluid;name)

    @named inport = CoolantPort()
    @named outport = CoolantPort()

    para = @parameters begin
        f = 0.03, [description = "Friction factor for laminar flow : f = 64/Re (Re < 2100)"]
        L , [description = "Length of the pipe"]
        D , [description = "Diameter of the pipe"]
    end

    vars = @variables begin
        P(t)
        s_in(t), [description = "Entropy of the fluid (J/kgK)"]
        p_in(t), [description = "Pressure of the fluid (Pa)"]
        T_in(t), [description = "Temperature of the fluid (K)"]
        h_in(t), [description = "Enthalpy of the fluid (J/kg)"]
        ρ_in(t), [description = "Density of the fluid (kg/m^3)"]

        s_out(t), [description = "Entropy of the fluid (J/kgK)"]
        p_out(t), [description = "Pressure of the fluid (Pa)"]
        T_out(t), [description = "Temperature of the fluid (K)"]
        h_out(t), [description = "Enthalpy of the fluid (J/kg)"]
        ρ_out(t), [description = "Density of the fluid (kg/m^3)"]

        Δp(t), [description = "Pressure drop across the pipe (Pa)"]
        mdot(t), [description = "Mass flow rate (kg/s)"]
        A(t), [description = "Cross-sectional area of the pipe (m^2)"]
    end

    eqs = [
        mdot ~ abs(inport.mdot)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)

        A ~ π*D^2/4
        Δp ~ (L*f*mdot^2)/(2*D*ρ_in*A^2)

        p_out ~ p_in - Δp
        T_out ~ T_in

        h_out ~ PropsSI("H","T",T_out,"P",p_out,fluid)
        s_out ~ PropsSI("S","H",h_out,"P",p_out,fluid)
        ρ_out ~ PropsSI("D","H",h_out,"P",p_out,fluid)
        outport.mdot ~ mdot
        outport.p ~ p_out
        outport.h ~ h_out
    ]

    compose(System(eqs, t, vars, para;name), inport, outport)
end


@component function PipeClapyeron(fluid::EoSModel = set_fluid;name)

    @named inport = CoolantPort()
    @named outport = CoolantPort()

    para = @parameters begin
        f = 0.03, [description = "Friction factor for laminar flow : f = 64/Re (Re < 2100)"]
        L , [description = "Length of the pipe"]
        D , [description = "Diameter of the pipe"]
    end

    vars = @variables begin
        s_in(t), [description = "Entropy of the fluid (J/K)"]
        p_in(t), [description = "Pressure of the fluid (Pa)"]
        T_in(t), [description = "Temperature of the fluid (K)"]
        h_in(t), [description = "Enthalpy of the fluid (J)"]
        ρ_in(t), [description = "Density of the fluid (kg/m^3)"]
        x_in(t), [description = "Mass fraction of first component in the mixture"]
        z_in(t), [description = "Mole fraction of the mixture"]

        s_out(t), [description = "Entropy of the fluid (J/K)"]
        p_out(t), [description = "Pressure of the fluid (Pa)"]
        T_out(t), [description = "Temperature of the fluid (K)"]
        h_out(t), [description = "Enthalpy of the fluid (J)"]
        ρ_out(t), [description = "Density of the fluid (kg/m^3)"]
        x_out(t), [description = "Mass fraction of first component in the mixture"]
        z_out(t), [description = "Mole fraction of the mixture"]

        Δp(t), [description = "Pressure drop across the pipe (Pa)"]
        mdot(t), [description = "Mass flow rate (g/s)"]
        A(t), [description = "Cross-sectional area of the pipe (m^2)"]
    end

    eqs = [
        x_in ~ inport.x
        z_in ~ mass_to_moles(fluid,x_in,mdot)
        mdot ~ abs(inport.mdot)
        p_in ~ inport.p
        T_in ~ ph_temperature(fluid,inport.p,inport.h,z_in)
        h_in ~ inport.h
        s_in ~ ph_entropy(fluid,inport.p,inport.h,z_in)
        ρ_in ~ ph_mass_density(fluid,inport.p,inport.h,z_in)

        A ~ π*D^2/4
        Δp ~ (L*f*(mdot/1000)^2)/(2*D*ρ_in*A^2)

        p_out ~ p_in - Δp
        T_out ~ T_in

        h_out ~ pt_enthalpy(fluid,p_out,T_out,z_in)
        s_out ~ ph_entropy(fluid,p_out,h_out,z_in)
        ρ_out ~ ph_mass_density(fluid,p_out,h_out,z_in)
        outport.mdot ~ inport.mdot
        outport.p ~ inport.p
        outport.h ~ inport.h
        outport.x ~ inport.x
    ]

    compose(System(eqs, t, vars, para;name), inport, outport)
end


"""
`Pipe(fluid::AbstractString = set_fluid;name)`

pressure drop across pipe using Darcy-Weisbach equation
"""
function Pipe(;name,fluid = set_fluid)
    if fluid isa EoSModel
        return PipeClapyeron(fluid;name = name)
    end
    if fluid isa AbstractString
        return PipeCoolProp(fluid;name = name)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
end


export Pipe
