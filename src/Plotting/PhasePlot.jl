function PhasePlot_TS(fig,fluid::EoSModel,p,z;N=100)
    p_min = minimum(p)
    p_max = CarnotCycles.CriticalPressure(fluid,z) 
    p_ = collect(range(p_min,p_max,N))
    T_bubble = zeros(N)
    T_dew = zeros(N)
    s_bubble = zeros(N)
    s_dew = zeros(N)
    for i in 1:N
        T_bubble[i] = CarnotCycles.qp_temperature(fluid,0.0,p_[i],z)
        s_bubble[i] = CarnotCycles.qp_entropy(fluid,0.0,p_[i],z)
        T_dew[i] = CarnotCycles.qp_temperature(fluid,1.0,p_[i],z)
        s_dew[i] = CarnotCycles.qp_entropy(fluid,1.0,p_[i],z)
    end
    plot!(fig,s_bubble,T_bubble,label = "Bubble Point",xlabel = "Entropy (J/K)", ylabel = "Temperature (K)",alpha = 0.3,linewidth=5)
    plot!(fig,s_dew,T_dew,label = "Dew Point",xlabel = "Entropy (J/K)", ylabel = "Temperature (K)",alpha = 0.3,linewidth=5)
end

function PhasePlot_TS(fig,fluid::AbstractString,p;N=100)
    p_min = minimum(p)
    p_max = CoolProp.PropsSI("PCRIT",fluid)
    p_ = collect(range(p_min,p_max,N))
    T_bubble = zeros(N)
    T_dew = zeros(N)
    s_bubble = zeros(N)
    s_dew = zeros(N)
    for i in 1:N
        T_bubble[i] = CarnotCycles.PropsSI("T","P",p_[i],"Q",0.0,fluid)
        s_bubble[i] = CarnotCycles.PropsSI("S","P",p_[i],"Q",0.0,fluid)
        T_dew[i] = CarnotCycles.PropsSI("T","P",p_[i],"Q",1.0,fluid)
        s_dew[i] = CarnotCycles.PropsSI("S","P",p_[i],"Q",1.0,fluid)
    end
    plot!(fig,s_bubble,T_bubble,label = "Bubble Point",xlabel = "Entropy (J/K)", ylabel = "Temperature (K)",linewidth=5,alpha = 0.3)
    plot!(fig,s_dew,T_dew,label = "Dew Point",xlabel = "Entropy (J/K)", ylabel = "Temperature (K)",linewidth=5,alpha = 0.3)
end

function PhasePlot_PH(fig,fluid::EoSModel,p,z,N=100)
    p_min = minimum(p)
    p_max = CarnotCycles.CriticalPressure(fluid,z) 
    p_ = collect(range(p_min,p_max,N))
    h_bubble = zeros(N)
    h_dew = zeros(N)
    for i in 1:N
        h_bubble[i] = CarnotCycles.qp_enthalpy(fluid,0.0,p_[i],z)
        h_dew[i] = CarnotCycles.qp_enthalpy(fluid,1.0,p_[i],z)
    end
    plot!(fig,h_bubble,p_,label = "Bubble Point",ylabel = "Pressure (Pa)", xlabel = "Enthalpy (J)",alpha = 0.3,linewidth=5)
    plot!(fig,h_dew,p_,label = "Dew Point",ylabel = "Pressure (Pa)", xlabel = "Enthalpy (J)",alpha = 0.3,linewidth=5)
    
end

function PhasePlot_PH(fig,fluid::AbstractString,p,N=100)
    p_min = minimum(p)
    p_max = CoolProp.PropsSI("PCRIT",fluid)
    p_ = collect(range(p_min,p_max,N))
    h_bubble = zeros(N)
    h_dew = zeros(N)
    for i in 1:N
        h_bubble[i] = CarnotCycles.PropsSI("H","P",p_[i],"Q",0.0,fluid)
        h_dew[i] = CarnotCycles.PropsSI("H","P",p_[i],"Q",1.0,fluid)
    end
    plot!(fig,h_bubble,p_,label = "Bubble Point",xlabel = "Enthalpy (J/kg)", ylabel = "Pressure (Pa)",linewidth=5,alpha = 0.3)
    plot!(fig,h_dew,p_,label = "Dew Point",xlabel = "Enthalpy (J/kg)", ylabel = "Pressure (Pa)",linewidth=5,alpha = 0.3)
    
end

function plot_PH(sol::SteadyStateSolution,system::Vector{ODESystem},names::Vector{String};phase = true,fluid = CarnotCycles.set_fluid)
    if fluid isa EoSModel
        @assert length(system)-2 == length(names) "Length of system and names must be equal, excluding source and sink"
        system_ = system[2:end-1] # Exclude source and sink
        z = sol[system_[1].z_in]
        h = zeros(length(system_))
        p = zeros(length(system_))
        for i in eachindex(names)
            h[i] = sol[system_[i].h_in]
            p[i] = sol[system_[i].p_in]
        end    
        append!(h,sol[system[end].h])
        append!(p,sol[system[end].p])
        fig = Plots.plot()
    
        if phase == true
            PhasePlot_PH(fig,fluid,p,z)
        end
        for i = 1: length(names)
            h_ = collect(range(h[i],h[i+1],100))
            p_ = collect(range(p[i],p[i+1],100))
            fig = plot!(fig,h_,p_,label = names[i],xlabel = "Enthalpy (J/K)", ylabel = "Pressure (Pa)", title = "Pressure-Enthalpy Diagram")
        end
        return fig
    end

    if fluid isa AbstractString
        @assert length(system)-2 == length(names) "Length of system and names must be equal, excluding source and sink"
        system_ = system[2:end-1] # Exclude source and sink
        h = zeros(length(system_))
        p = zeros(length(system_))
        for i in eachindex(names)
            h[i] = sol[system_[i].h_in]
            p[i] = sol[system_[i].p_in]
        end    
        append!(h,sol[system[end].h])
        append!(p,sol[system[end].p])
        fig = Plots.plot()
    
        if phase == true
            PhasePlot_PH(fig,fluid,p)
        end
        for i = 1: length(names)
            h_ = collect(range(h[i],h[i+1],100))
            p_ = collect(range(p[i],p[i+1],100))
            fig = plot!(fig,h_,p_,label = names[i],xlabel = "Enthalpy (J/K/kg)", ylabel = "Pressure (Pa)", title = "Pressure-Enthalpy Diagram")
        end
        return fig
    end
end

function plot_TS(sol::SteadyStateSolution,system::Vector{ODESystem},names::Vector{String};phase = true,fluid = CarnotCycles.set_fluid)
    if fluid isa EoSModel
        @assert length(system)-2 == length(names) "Length of system and names must be equal, excluding source and sink"
        system_ = system[2:end-1] # Exclude source and sink
        z = sol[system_[1].z_in]
        h = zeros(length(system_))
        p = zeros(length(system_))
        for i in eachindex(names)
            h[i] = sol[system_[i].h_in]
            p[i] = sol[system_[i].p_in]
        end    
        append!(h,sol[system[end].h])
        append!(p,sol[system[end].p])
        fig = Plots.plot()
    
        if phase == true
            PhasePlot_TS(fig,fluid,p,z)
        end
        for i = 1: length(names)
            h_ = collect(range(h[i],h[i+1],100))
            p_ = collect(range(p[i],p[i+1],100))
            T_ = similar(h_)
            s_ = similar(h_)
            for j in eachindex(h_)
                T_[j] = CarnotCycles.ph_temperature_qp(fluid,p_[j],h_[j],CarnotCycles.mass_to_moles(CarnotCycles.set_fluid,sol[system_[i].inport.x],sol[system_[i].inport.mdot]))
                s_[j] = CarnotCycles.ph_entropy_qp(fluid,p_[j],h_[j],CarnotCycles.mass_to_moles(CarnotCycles.set_fluid,sol[system_[i].inport.x],sol[system_[i].inport.mdot]))
            end
            fig = plot!(fig,s_,T_,label = names[i],xlabel = "Entropy (J/K)", ylabel = "Temperature (K)", title = "Temperature-Entropy Diagram")
        end
        return fig
    end

    if fluid isa AbstractString
        @assert length(system)-2 == length(names) "Length of system and names must be equal, excluding source and sink"
        system_ = system[2:end-1] # Exclude source and sink
        h = zeros(length(system_))
        p = zeros(length(system_))
        for i in eachindex(names)
            h[i] = sol[system_[i].h_in]
            p[i] = sol[system_[i].p_in]
        end    
        append!(h,sol[system[end].h])
        append!(p,sol[system[end].p])
        fig = Plots.plot()
    
        if phase == true
            PhasePlot_TS(fig,fluid,p)
        end
        for i = 1: length(names)
            h_ = collect(range(h[i],h[i+1],100))
            p_ = collect(range(p[i],p[i+1],100))
            T_ = similar(h_)
            s_ = similar(h_)
            for j in eachindex(h_)
                T_[j] = CarnotCycles.PropsSI("T","P",p_[j],"H",h_[j],fluid)
                s_[j] = CarnotCycles.PropsSI("S","P",p_[j],"H",h_[j],fluid)
            end
            fig = plot!(fig,s_,T_,label = names[i],xlabel = "Entropy (J/K/kg)", ylabel = "Temperature (K)", title = "Temperature-Entropy Diagram")
        end
        return fig
    end
end

"""
`plot(sol::SteadyStateSolution,system::Vector{ODESystem},names::Vector{String};phase = true,fluid = CarnotCycles.set_fluid,type = :TS)`
    * Plots the phase diagram of the system using the given solution and system.
    * `sol`: The solution object containing the results of the simulation. Do not include source and sink in the system.
    * `system`: The vector of ODESystem objects representing the system.
    * `names`: The vector of names for each component in the system.
    * `phase`: A boolean indicating whether to plot the phase boundaries or not. Default is true.
    * `fluid`: The fluid model to be used for plotting. Default is `CarnotCycles.set_fluid`.
    * `type`: The type of plot to be generated. Can be either :TS or :PH. Default is :TS.

    Returns a `Plots.plot` object.
"""
function plot(sol::SteadyStateSolution,system::Vector{ODESystem},names::Vector{String};phase = true,fluid = CarnotCycles.set_fluid,type = :TS)
    if type == :TS
        return plot_TS(sol,system,names;phase = phase,fluid = fluid)
    elseif type == :PH
        return plot_PH(sol,system,names;phase = phase,fluid = fluid)
    else
        error("Type must be either :TS or :PH")
    end
    
end

export PhasePlot_TS, PhasePlot_PH, plot_TS, plot_PH, plot