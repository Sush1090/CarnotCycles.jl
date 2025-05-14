
"""
`Pipe(fluid::AbstractString = set_fluid;name)`

pressure drop across pipe using Darcy-Weisbach equation
"""
@component function Pipe(fluid::AbstractString = set_fluid;name)

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
        h_out ~ h_in

        T_out ~ PropSI("T","H",h_out,"P",p_out,fluid)
        s_out ~ PropsSI("S","H",h_out,"P",p_out,fluid)
        ρ_out ~ PropsSI("D","H",h_out,"P",p_out,fluid)
        outport.mdot ~ mdot
        outport.p ~ p_out
        outport.h ~ h_out
    ]

    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end


