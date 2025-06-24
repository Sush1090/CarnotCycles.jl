"""
`IsentropicCompressor(;name,fluid=set_fluid)`

A compressor with isentropic Effeciency and pressure ratio as a parameter is chosen. 
"""
@component function IsentropicCompressor(;name,fluid=set_fluid)
    if fluid isa AbstractString
        return IsentropicCompressorCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return IsentropicCompressorClapeyron(name=name,fluid=fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
end


@component function IsentropicCompressorCoolProp(;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)
    @named power =  PowerPort()
    para = @parameters begin
        η=0.8, [description = "Isentropic Effeciency",bounds = (0.0,1.0)]
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
            outport.p ~ πc * inport.p
            outport.h ~ IsentropicCompression(πc, inport.h, inport.p,fluid,η)
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
            power.P ~ P
   ]
   compose(System(eqs, t, vars, para;name), inport, outport,power)
end


function IsentropicCompressorClapeyron(;name,fluid = set_fluid) 
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)
    @named power =  PowerPort()
    para = @parameters begin
        η=0.8, [description = "Isentropic Effeciency",bounds = (0.0,1.0)]
        πc, [description = "Pressure ratio"]
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
            z_in ~ mass_to_moles(fluid,x_in,mdot_in)
            p_in ~ inport.p
            h_in ~ inport.h
            s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
            T_in ~ ph_temperature(fluid,p_in,h_in,z_in)
            ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)

            h_out ~ IsentropicCompressionClapeyron(πc,h_in,p_in,z_in,fluid,η)
            p_out ~ πc*p_in
            s_out ~ ph_entropy(fluid,p_out,h_out,z_out)
            z_out ~ z_in
            T_out ~ ph_temperature(fluid,p_out,h_out,z_out)
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
   compose(System(eqs, t, vars, para;name), inport, outport,power)
end

export IsentropicCompressor



function IsothermalCompressorCoolProp(;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    # @unpack πc = type
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)
    @named power =  PowerPort()
    para = @parameters begin
        πc, [description = "Pressure Ratio"]
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
            outport.p ~ πc * inport.p
            outport.h ~ IsothermalCompression(πc, inport.h, inport.p,fluid)
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
            power.P ~ P
   ]
   compose(System(eqs, t, vars, para;name), inport, outport,power)
end

function IsothermalCompressorClapeyron(;name,fluid = set_fluid) 
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)
    @named power =  PowerPort()
    para = @parameters begin
        πc, [description = "Pressure Ratio"]
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)
        z_in(t)
        x_in(t)
        mdot_in(t)

        mdot_out(t)
        x_out(t)
        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
        z_out(t)
     end
   eqs = [  mdot_in ~ inport.mdot
            x_in ~ inport.x
            z_in ~ mass_to_moles(fluid,x_in,mdot_in)
            p_in ~ inport.p
            h_in ~ inport.h
            s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
            ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)
            T_in ~ ph_temperature(fluid,p_in,h_in,z_in)

            h_out ~ IsothermalCompressionClapeyron(πc,h_in,p_in,z_in,fluid)
            p_out ~ p_in*πc
            s_out ~ ph_entropy(fluid,p_out,h_out,z_out)
            T_out ~ ph_temperature(fluid,p_out,h_out,z_out)
            ρ_out ~ ph_mass_density(fluid,p_out,h_out,z_out)
            z_out ~ z_in
            x_out ~ x_in
            mdot_out ~ mdot_in

            P ~ h_out - h_in

            outport.p ~ p_out
            outport.h ~ h_out
            outport.x ~ x_out
            outport.mdot ~ mdot_out
            power.P ~ P
   ]
   compose(System(eqs, t, vars, para;name), inport, outport,power)
end

"""
`IsothermalCompressor(;name,fluid = set_fluid)`

A compressor with pressure ratio as a parameter is chosen. 
"""
@component function IsothermalCompressor(;name,fluid = set_fluid)
    if fluid isa AbstractString
        return IsothermalCompressorCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return IsothermalCompressorClapeyron(name=name,fluid=fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected")) 
    end
end
export IsothermalCompressor




function IsochoricCompressorCoolProp(;name,fluid = set_fluid) 
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)
    @named power =  PowerPort()
    para = @parameters begin
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
            outport.p ~ πc * inport.p
            outport.h ~ IsochoricCompression(πc, inport.h, inport.p,fluid)
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
            power.P ~ P 
   ]
   compose(System(eqs, t, vars, para;name), inport, outport,power)
end


function IsochoricCompressorClapeyron(;name,fluid = set_fluid) 
    @named inport = CoolantPort(fluid=fluid)
    @named outport = CoolantPort(fluid=fluid)
    @named power =  PowerPort()
    para = @parameters begin
        πc, [description = "Pressure ratio"]
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)
        x_in(t)
        mdot_in(t)
        z_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
        x_out(t)
        mdot_out(t)
        z_out(t)
     end
   eqs = [  mdot_in ~ inport.mdot
            x_in ~ inport.x
            z_in ~ mass_to_moles(fluid,x_in,mdot_in)
            p_in ~ inport.p
            h_in ~ inport.h
            s_in ~ ph_entropy(fluid,p_in,h_in,z_in)
            ρ_in ~ ph_mass_density(fluid,p_in,h_in,z_in)
            T_in ~ ph_temperature(fluid,p_in,h_in,z_in)

            h_out ~ IsochoricCompressionClapeyron(πc,h_in,p_in,z_in,fluid)
            p_out ~ p_in*πc
            s_out ~ ph_entropy(fluid,p_out,h_out,z_out)
            T_out ~ ph_temperature(fluid,p_out,h_out,z_out)
            ρ_out ~ ph_mass_density(fluid,p_out,h_out,z_out)
            z_out ~ z_in
            x_out ~ x_in
            mdot_out ~ mdot_in

            P ~ h_out - h_in

            outport.p ~ p_out
            outport.h ~ h_out
            outport.x ~ x_out
            outport.mdot ~ mdot_out

            power.P ~ P
   ]
   compose(System(eqs, t, vars, para;name), inport, outport,power)
end

"""
`IsochoricCompressor(;name,fluid = set_fluid)`

A compressor with pressure ratio as a parameter is chosen. 
"""
@component function IsochoricCompressor(;name,fluid = set_fluid)
    if fluid isa AbstractString
        return IsochoricCompressorCoolProp(name=name,fluid=fluid)
    end
    if fluid isa EoSModel
        return IsochoricCompressorClapeyron(name=name,fluid=fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected")) 
    end
end

export IsochoricCompressor


