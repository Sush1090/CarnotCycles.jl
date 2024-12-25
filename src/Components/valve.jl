
"""
`Valve(;name,fluid): Isenthalpic valve`
"""
function Valve(;name,πc,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    if fluid isa AbstractString
        return ValveCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return ValveClpeyron(name=name,fluid=fluid)
    end
end


function ValveCoolProp(;name,fluid=set_fluid)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
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
    para = @parameters begin
        πc,[description = "Pressure ratio (-)"]
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot) 
        πc ~  inport.p/outport.p
        outport.h ~ inport.h
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


@component function ValveClapeyron(;name,fluid=set_fluid)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)
        mdot_in(t)
        x_in(t)
        z_in(t)

        z_out(t)
        mdot_out(t)
        x_out(t)
        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
    end
    para = @parameters begin
        πc,[description = "Pressure ratio (-)"]
    end
    eqs = [
        mdot_in ~ inport.mdot
        z_in ~ mass_to_moles(fluid,[x_in,1-x_in],mdot_in)
        p_in ~ inport.p
        h_in ~ inport.h
        s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
        T_in ~ ph_temperature(fluid,p_in,h_in,z_in)
        ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)
        x_in ~ inport.x

        mdot_out ~ mdot_in
        z_out ~ z_in
        p_out ~ p_in/πc
        h_out ~ h_in
        s_out ~ ph_entropy(fluid,p_out,h_out,z_out)
        T_out ~ ph_temperature(fluid,p_out,h_out,z_out)
        ρ_out ~ ph_mass_density(fluid,p_out,h_out,z_out)
        x_out ~ x_in

        outport.h ~ h_out
        outport.p ~ p_out
        outport.x ~ x_out
        outport.mdot ~ mdot_out

    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end


# struct ThreeFacedValveSplit
# ratio::AbstractVector
# function ThreeFacedValveSplit(;ratio::AbstractVector)
#     @assert isapprox(sum(ratio),1) "Sum of the ratio of mass flow rates should be 1"
#     @assert size(ratio,1) == 2 "Only splits in two streams hence provide only two ratios"
#     new(ratio)
# end
# end


# """
# `Valve(type::ThreePhaseValveSplit;name,fluid=set_fluid)`

# """
# function Valve(type::ThreeFacedValveSplit;name,fluid=set_fluid)
#     if isnothing(fluid)
#         throw(error("Fluid not selected"))
#     end
#     @unpack ratio= type
#     @named inport = CoolantPort()
#     @named outport1 = CoolantPort()
#     @named outport2 = CoolantPort()
#     vars = @variables begin
#         P1(t)
#         P2(t)
#         P(t)

#         s_in(t)
#         p_in(t)
#         T_in(t)
#         h_in(t)
#         ρ_in(t)

#         s_out1(t)
#         p_out1(t)
#         T_out1(t)
#         h_out1(t)
#         ρ_out1(t)

#         s_out2(t)
#         p_out2(t)
#         T_out2(t)
#         h_out2(t)
#         ρ_out2(t)
#     end
#     para = @parameters begin
        
#     end
#     eqs = [
#         outport1.mdot ~ abs(inport.mdot)*ratio[1]
#         outport2.mdot ~ abs(inport.mdot)*ratio[2]
#         outport1.p ~  inport.p
#         outport2.p ~  inport.p
#         outport1.h ~ inport.h
#         outport2.h ~ inport.h
#         P1 ~ abs(outport1.mdot)*(outport1.h - inport.h)
#         P2 ~ abs(outport2.mdot)*(outport2.h - inport.h)
#         P ~ P1+ P2
#         s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
#         p_in ~ inport.p
#         T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
#         h_in ~ inport.h
#         ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)

#         s_out1 ~ PropsSI("S","H",outport1.h,"P",outport1.p,fluid)
#         p_out1 ~ outport1.p
#         T_out1 ~ PropsSI("T","H",outport1.h,"P",outport1.p,fluid)
#         h_out1 ~ outport1.h
#         ρ_out1 ~ PropsSI("D","H",outport1.h,"P",outport1.p,fluid)

#         s_out2 ~ PropsSI("S","H",outport2.h,"P",outport2.p,fluid)
#         p_out2 ~ outport2.p
#         T_out2 ~ PropsSI("T","H",outport2.h,"P",outport2.p,fluid)
#         h_out2 ~ outport2.h
#         ρ_out2 ~ PropsSI("D","H",outport2.h,"P",outport2.p,fluid)
#     ]
#     compose(ODESystem(eqs, t, vars, para;name), inport, outport1,outport2)
# end



# struct ThreeFacedValveCombine

# end


# """
# `Valve(type::ThreePhaseValveSplit;name,fluid=set_fluid)`

# """
# function Valve(type::ThreeFacedValveCombine;name,fluid=set_fluid)
#     if isnothing(fluid)
#         throw(error("Fluid not selected"))
#     end
#     @named inport1 = CoolantPort()
#     @named inport2 = CoolantPort()
#     @named outport = CoolantPort()
#     vars = @variables begin
#         r(t)

#         P1(t)
#         P2(t)
#         P(t)

#         s_in2(t)
#         p_in2(t)
#         T_in2(t)
#         h_in2(t)
#         ρ_in2(t)

#         s_in1(t)
#         p_in1(t)
#         T_in1(t)
#         h_in1(t)
#         ρ_in1(t)

#         s_out(t)
#         p_out(t)
#         T_out(t)
#         h_out(t)
#         ρ_out(t)
#     end
#     para = @parameters begin
        
#     end
#     eqs = [
#         outport.mdot ~ abs(inport1.mdot) + abs(inport2.mdot) 
#         outport.p ~  inport1.p*(abs(inport1.mdot)/outport.mdot) + inport2.p*(abs(inport2.mdot)/outport.mdot) 
#         outport.h ~ inport1.h*(abs(inport1.mdot)/outport.mdot) + inport2.h*(abs(inport2.mdot)/outport.mdot)
        
#         P1 ~ (abs(inport1.mdot)/outport.mdot)*(outport.h - inport1.h)
#         P2 ~ (abs(inport2.mdot)/outport.mdot)*(outport.h - inport2.h)
#         P ~ P1+ P2
        
#         s_in1 ~ PropsSI("S","H",inport1.h,"P",inport1.p,fluid)
#         p_in1 ~ inport1.p
#         T_in1 ~ PropsSI("T","H",inport1.h,"P",inport1.p,fluid)
#         h_in1 ~ inport1.h
#         ρ_in1 ~ PropsSI("D","H",inport1.h,"P",inport1.p,fluid)

#         s_in2 ~ PropsSI("S","H",inport2.h,"P",inport2.p,fluid)
#         p_in2 ~ inport2.p
#         T_in2 ~ PropsSI("T","H",inport2.h,"P",inport2.p,fluid)
#         h_in2 ~ inport2.h
#         ρ_in2 ~ PropsSI("D","H",inport2.h,"P",inport2.p,fluid)


#         s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
#         p_out ~ outport.p
#         T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
#         h_out ~ outport.h
#         ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
#     ]
#     compose(ODESystem(eqs, t, vars, para;name), inport1, inport2,outport)
# end



export Valve, ValveClapeyron, ValveCoolProp