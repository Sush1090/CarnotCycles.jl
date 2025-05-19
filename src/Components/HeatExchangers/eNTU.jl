"""
`flow` can be either `:counterflow`, `:parallelflow`, or `:crossflow`

returns ϵ for the given NTU and C_r and flow type   

C_r is the heat capacity ratio, defined as C_r = C_min / C_max, where C_min and C_max are the minimum and maximum heat capacities of the two fluids in the heat exchanger multiplied with their mass flow rate.
"""
function eNTU(NTU,C_r,flow::Symbol = :counterflow)
    if flow == :counterflow
        return (1 - exp(-NTU * (1 - C_r))) / (1 - C_r*exp(-NTU * (1 - C_r)))   
    end
    if flow == :parallelflow
        return (1 - exp(-NTU * (1 + C_r))) / (1 + C_r) 
    end

    if flow == :crossflow
        return 1 - exp((1 / C_r)*(NTU^0.22)*(exp(-C_r*(NTU^0.78)) - 1))
    end

    throw(ArgumentError("flow must be :counterflow, :parallelflow, or :crossflow. Check the flow type"))
end


function C_min(mdot_wf::Float64,cp_wf::Float64,mdot_sf::Float64,cp_sf::Float64)
    C_min = min(mdot_wf*cp_wf,mdot_sf*cp_sf)
    return C_min
end


function C_r(mdot_wf::Float64,cp_wf::Float64,mdot_sf::Float64,cp_sf::Float64)
    C_min = min(mdot_wf*cp_wf,mdot_sf*cp_sf)
    C_max = max(mdot_wf*cp_wf,mdot_sf*cp_sf)
    C_r = C_min / C_max
end

function NTU(hex::ϵNTU,C_min)
    @unpack flow, U, A = hex
    return (U*A)/C_min
end

function compute_Qdot_max(mdot_in_wf,T_in_wf,cp_in_wf,mdot_in_htf,T_in_htf,cp_in_htf)
    Qdot_max = C_min(mdot_in_wf,cp_in_wf,mdot_in_htf,cp_in_htf) * (T_in_wf - T_in_htf)
    return Qdot_max
end

function compute_Qdot(hex::HEXModel,mdot_in_1,T_in_1,cp_in_1,mdot_in_2,T_in_2,cp_in_2)
    Qdot_max = compute_Qdot_max(mdot_in_1,T_in_1,cp_in_1,mdot_in_2,T_in_2,cp_in_2)
    C_r = C_r(mdot_in_1,cp_in_1,mdot_in_2,cp_in_2)
    C_min = C_min(mdot_in_1,cp_in_1,mdot_in_2,cp_in_2)
    NTU = NTU(hex,C_min)
    ϵ = eNTU(NTU,C_r,hex.flow)
    Qdot = ϵ * Qdot_max
    return Qdot
end

@register_symbolic compute_Qdot(hex::HEXModel,mdot_in_1,T_in_1,cp_in_1,mdot_in_2,T_in_2,cp_in_2)


function HeatExchanger(model::HEXModel;name,fluid = set_fluid)
    if fluid isa EoSModel
        return HeatExchangerClapeyron(model;name = name,fluid = fluid)
    end
    if fluid isa AbstractString
        return HeatExchangerCoolProp(model;name = name,fluid = fluid)
    end
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
end

function HeatExchangerCoolProp(hex::ϵNTU;name = name,fluid)
    @named inport = CoolantPort()
    @named outport = CoolantPort()

    @named htf_inport = StoragePort()
    @named htf_outport = StoragePort()

    para = @parameters begin
        Cp_htf(t), [description = "Specific heat of HTF (J/kg-K)"]
        T_htf_in(t), [description = "Inlet Temperature of HTF"]
        mdot_htf(t), [description = "Inlet mass flow rate of HTF"]
    end

    vars = @variables begin
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)
        mdot_in(t)

        Cp_in(t), [description = "Specific heat of HTF (J/kg-K)"]

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
        mdot_out(t)
        cp_out(t), [description = "Specific heat of HTF (J/kg-K)"]

        Qdot(t), [description = "Heat transfer rate (W)"]
        ϵ(t), [description = "Effectiveness"]

        T_htf_out(t), [description = "Outlet Temperature of HTF"]
    end

    eqs = [
        h_in ~ inport.h
        p_in ~ inport.p
        mdot_in ~ inport.mdot
        T_in ~ PropsSI("T","H",h_in,"P",p_in,fluid)
        s_in ~ PropsSI("S","H",h_in,"P",p_in,fluid)
        ρ_in ~ PropsSI("D","H",h_in,"P",p_in,fluid)
        Cp_in ~ PropsSI("CPMASS","H",h_in,"P",p_in,fluid)

        htf_inport.mdot ~ mdot_htf
        htf_inport.T ~ T_htf_in
        
        Qdot ~ compute_Qdot(hex,mdot_in,T_in,Cp_in,mdot_htf,T_htf_in,Cp_htf)
        p_out ~ p_in
        mdot_out ~ mdot_in

        Qdot ~ outport.mdot * (h_out - h_in)
        T_out ~ PropsSI("T","H",h_out,"P",p_out,fluid)
        s_out ~ PropsSI("S","H",h_out,"P",p_out,fluid)
        ρ_out ~ PropsSI("D","H",h_out,"P",p_out,fluid)
        Cp_out ~ PropsSI("CPMASS","H",h_out,"P",p_out,fluid)
        p_out ~ outport.p
        outport.mdot ~ mdot_out

        T_htf_out ~ T_htf_in - Qdot/(mdot_htf * Cp_htf)
        htf_outport.T ~ T_htf_out
        htf_outport.mdot ~ mdot_htf
    ]

    compose(ODESystem(eqs, t, vars, para;name), inport, outport, htf_inport, htf_outport)
end