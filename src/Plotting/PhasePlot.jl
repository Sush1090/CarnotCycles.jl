

abstract type PhasePlotType end

struct PhasePlotType_PH <: PhasePlotType end
struct PhasePlotType_TS <: PhasePlotType end

export PhasePlotType_PH,PhasePlotType_TS

T_line(p,h,fluid) = PropsSI("T","P",p,"H",h,fluid)
S_line(p,h,fluid) = PropsSI("S","P",p,"H",h,fluid)
ρ_line(p,h,fluid) = PropsSI("D","P",p,"H",h,fluid)

function LineData(f::Function,x1,y1,x2,y2,fluid::AbstractString=set_fluid;points = 100)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    x = collect(range(x1,x2,points))
    y = collect(range(y1,y2,points))
    data = f.(x,y,fluid)
    return x,y,data
end


function CollectPhaseData(sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid;points = 100,p_width::Tuple = (1e3,1e3), h_width::Tuple = (1e3,1e3))
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    p_crit = CoolProp.PropsSI("PCRIT",fluid)
    plot_sys = system[2:end-1];
    propx = :p
    propy = :h
    pp = [sol[getproperty(i.inport, propx)][1] for i in plot_sys]
    hh = [sol[getproperty(i.inport, propy)][1] for i in plot_sys]
    pmin = minimum(pp); pmax = maximum(pp)
    hmin = minimum(hh); hmax = maximum(hh)
    p_range = collect(range(pmin-p_width[1],pmax + p_width[2],points));
    h_range = collect(range(hmin-h_width[1],hmax + h_width[2],points));

    p_sat_curve = collect(range(pmin,p_crit,points))
    T_sat = PropsSI.("T","Q",1,"P",p_sat_curve,fluid)
    s_sat_liquid = PropsSI.("S","Q",0,"P",p_sat_curve,fluid)
    s_sat_gas = PropsSI.("S","Q",1,"P",p_sat_curve,fluid)
    h_sat_liquid = PropsSI.("H","Q",0,"P",p_sat_curve,fluid)
    h_sat_gas = PropsSI.("H","Q",1,"P",p_sat_curve,fluid)
    return p_range,h_range,T_sat,s_sat_liquid,s_sat_gas,h_sat_liquid,h_sat_gas,p_sat_curve
end


function PhasePlot(type::PhasePlotType_TS,sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    plot_sys = system[2:end-1];
    propx = :p
    propy = :h
    pp = [sol[getproperty(i.inport, propx)][1] for i in plot_sys]
    hh = [sol[getproperty(i.inport, propy)][1] for i in plot_sys]

    p_range,h_range,T_sat,s_sat_liquid,s_sat_gas,h_sat_liquid,h_sat_gas,p_sat_curve = CollectPhaseData(sol,system,fluid)
    T_cycle_points = PropsSI.("T","P",pp,"H",hh,fluid) 
    s_cycle_points = PropsSI.("S","P",pp,"H",hh,fluid)
    T_data = []; s_data = []
    cycle_data_pp = [pp;pp[1]]; cycle_data_hh = [hh;hh[1]]
    for i = 1:size(pp,1)
        x1 = cycle_data_pp[i]; y1 = cycle_data_hh[i];
        x2 = cycle_data_pp[i+1]; y2 = cycle_data_hh[i+1];
        _,_,TT = LineData(T_line,x1,y1,x2,y2,fluid)
        T_data = append!(T_data,TT)
        _,_,ss = LineData(S_line,x1,y1,x2,y2,fluid)
        s_data = append!(s_data,ss)
    end
    fig1 = plot(s_data,T_data,label = "Cycle")
     scatter!(s_cycle_points,T_cycle_points,label="Cycle Points")
     plot!(s_sat_liquid,T_sat,label ="Saturation Curve liquid")
     plot!(s_sat_gas,T_sat,label ="Saturation Curve gas")
    xlabel!("Entropy (J/kg/K)")
    ylabel!("Temperature (K)")
    return fig1
end


# """
# `PhasePlot(type::PhasePlotType_PH,sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid)`

# Plots the saturation curve of the fluid along with the cycle state-points. 
#  * Arguments:
#     1. `type`   : Chose between `PhasePlotType_PH` or `PhasePlotType_TS`
#     2. `sol`    : The solution from the ODE system. 
#     3. `system` : The vector of components chosen. The first should be `MassSource` and the last has to be `MassSink`.
# """
function PhasePlot(type::PhasePlotType_PH,sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    plot_sys = system[2:end-1];
    propx = :p
    propy = :h
    pp = [sol[getproperty(i.inport, propx)][1] for i in plot_sys]
    hh = [sol[getproperty(i.inport, propy)][1] for i in plot_sys]
    p_range,h_range,T_sat,s_sat_liquid,s_sat_gas,h_sat_liquid,h_sat_gas,p_sat_curve = CollectPhaseData(sol,system,fluid)
    p_data = []; h_data = []
    cycle_data_pp = [pp;pp[1]]; cycle_data_hh = [hh;hh[1]]
    for i = 1:size(pp,1)
        x1 = cycle_data_pp[i]; y1 = cycle_data_hh[i];
        x2 = cycle_data_pp[i+1]; y2 = cycle_data_hh[i+1];
        h_ = collect(range(y1,y2,100)); h_data = append!(h_data,h_)
        p_ = collect(range(x1,x2,100)); p_data = append!(p_data,p_)
    end
    fig1 = plot(h_sat_liquid,p_sat_curve,label = "Saturation Curve liquid")
    plot!(h_sat_gas,p_sat_curve,label = "Saturation Curve gas")
    scatter!(hh,pp,label="Cycle Points")
    xlabel!("Enthalpy (J/kg)")
    ylabel!("Pressure (Pa)")
    plot!(h_data,p_data,label = "Cycle")
    return fig1
end

export PhasePlot

# """
# `CyclePlot(type::PhasePlotType_PH,sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid)`

# Added Cycle ploting without saturation curve. Only plots the  cycle state points and their path on P-H plane. 

#  * Arguments:
#     1. `type`   : Chose between `PhasePlotType_PH` or `PhasePlotType_TS`
#     2. `sol`    : The solution from the ODE system. 
#     3. `system` : The vector of components chosen. The first should be `MassSource` and the last has to be `MassSink`.
#     4. `fluid`  : The fluid of chosen. Defaults to `set_fluid`
# """
function CyclePlot(type::PhasePlotType_PH,sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    plot_sys = system[2:end-1];
    propx = :p
    propy = :h
    pp = [sol[getproperty(i.inport, propx)][1] for i in plot_sys]
    hh = [sol[getproperty(i.inport, propy)][1] for i in plot_sys]
    p_data = []; h_data = []
    cycle_data_pp = [pp;pp[1]]; cycle_data_hh = [hh;hh[1]]
    for i = 1:size(pp,1)
        x1 = cycle_data_pp[i]; y1 = cycle_data_hh[i];
        x2 = cycle_data_pp[i+1]; y2 = cycle_data_hh[i+1];
        h_ = collect(range(y1,y2,100)); h_data = append!(h_data,h_)
        p_ = collect(range(x1,x2,100)); p_data = append!(p_data,p_)
    end
    fig1 =  scatter(hh,pp,label="Cycle Points")
    plot!(h_data,p_data,label = "Cycle")
    xlabel!("Enthalpy (J/kg)")
    ylabel!("Pressure (Pa)")
    return fig1
end

# """
# `CyclePlot(type::PhasePlotType_TS,sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid)`

# Added Cycle ploting without saturation curve. Only plots the  cycle state points and their path on T-S plane. 

#  * Arguments:
#     1. `type`   : Chose between `PhasePlotType_PH` or `PhasePlotType_TS`
#     2. `sol`    : The solution from the ODE system. 
#     3. `system` : The vector of components chosen. The first should be `MassSource` and the last has to be `MassSink`.
#     4. `fluid`  : The fluid of chosen. Defaults to `set_fluid`
# """
function CyclePlot(type::PhasePlotType_TS,sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString=set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    plot_sys = system[2:end-1];
    propx = :p
    propy = :h
    pp = [sol[getproperty(i.inport, propx)][1] for i in plot_sys]
    hh = [sol[getproperty(i.inport, propy)][1] for i in plot_sys]
    T_cycle_points = PropsSI.("T","P",pp,"H",hh,fluid) 
    s_cycle_points = PropsSI.("S","P",pp,"H",hh,fluid)
    T_data = []; s_data = []
    cycle_data_pp = [pp;pp[1]]; cycle_data_hh = [hh;hh[1]]
    for i = 1:size(pp,1)
        x1 = cycle_data_pp[i]; y1 = cycle_data_hh[i];
        x2 = cycle_data_pp[i+1]; y2 = cycle_data_hh[i+1];
        _,_,TT = LineData(T_line,x1,y1,x2,y2,fluid)
        T_data = append!(T_data,TT)
        _,_,ss = LineData(S_line,x1,y1,x2,y2,fluid)
        s_data = append!(s_data,ss)
    end
    fig1 = plot(s_data,T_data,label = "Cycle")
     scatter!(s_cycle_points,T_cycle_points,label="Cycle Points")
     xlabel!("Entropy (J/kg/K)")
     ylabel!("Temperature (K)")
     return fig1
end

export CyclePlot