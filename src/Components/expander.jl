
@component function IsentropicExpander(;name,fluid=set_fluid)
    if fluid isa AbstractString
        return IsentropicExpanderCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return IsentropicExpanderClapeyron(name=name,fluid = fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
end
export IsentropicExpander

function IsentropicExpanderCoolProp(;name,fluid=set_fluid)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        η, [description = "Isentropic Effeciency"]
        πc, [description = "Pressure ratio"]
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~  inport.p/πc
            outport.h ~ IsentropicExpansion(πc, inport.h, inport.p,fluid,η)
            P ~ abs(inport.mdot)*(outport.h - inport.h)
            s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
            p_in ~ inport.p
            T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
            h_in ~ inport.h
            ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
            p_out ~ outport.p
            T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
            h_out ~ outport.h
            ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

function IsentropicExpanderClapeyron(;name,fluid = set_fluid) 
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    @named power =  PowerPort()
    para = @parameters begin
        η=0.8, [description = "Isentropic Effeciency",bounds = (0.0,1.0)]
        πc, [description = "Pressure ratio"]
        pressure_corrector = 1 , [description = "Small pressure added for numerical stability"]
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)
        z_in(t)
        mdot_in(t)
        x_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
        z_out(t)
        mdot_out(t)
        x_out(t)
     end
   eqs = [  
            mdot_in ~ inport.mdot
            x_in ~ inport.x
            z_in ~ mass_to_moles(fluid,[x_in,1-x_in],mdot_in)
            p_in ~ inport.p
            h_in ~ inport.h
            s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
            T_in ~ ph_temperature(fluid,p_in,h_in,z_in)
            ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)

            h_out ~ pt_enthalpy(fluid,p_out,T_out,z_out)
            p_out ~ p_in/πc
            s_out ~ ph_entropy(fluid,p_out,h_out,z_out)
            z_out ~ z_in
            T_out ~ PT_IsentropicExpansionClapeyron(fluid,T_in,p_in+pressure_corrector,z_in,πc,η)
            ρ_out ~ ph_mass_density(fluid,p_out,h_out,z_out)
            mdot_out ~ mdot_in
            x_out ~ x_in   

            P ~ h_out - h_in

            outport.p ~ p_out
            outport.h ~ h_out
            outport.x ~ x_out
            outport.mdot ~ mdot_out

            power.P ~ P
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport,power)
end


