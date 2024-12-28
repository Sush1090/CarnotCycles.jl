
"""

"""
@component function PackedBed(;name,N)
    @named inlet = StoragePort()
    @named outlet = StoragePort()
    vars = @variables begin
        Ts(t)[1:N+1], [description = "Solid Phase Temperature (K)"]
        Tg(t)[1:N+1], [description = "Gas Phase Temperature (K)"]
        x(t)[1:N+1], [description = "Length along the axis (m)"]
        Δx(t), [description = "Descritized length (m)"]
        Vdot(t), [description = "Volumetric Flow Rate (m³/s)"]
    end
    para = @parameters begin
        ϵ = 0.42, [description = "Void Fraction (-)"]
        ρs = 2801, [description = "Density of solid (kg/m³)"]
        ρg = 869, [description = "Density of gas (kg/m³)"]
        T_amb=298, [description = "Ambient Temperature (K)"]
        Cg = 1800, [description = "Specific Capacity of Fluid (J/kg/K)"]
        Cs = 745, [description = "Specific Capacity of Solid (J/kg/K)"]
        A = 0.067, [description = "Cross Section Area of storage (m²)"]
        UA = 5, [description = "Heat loss coeffecient (W/K)"]
        # V =0.02, [description = "Volume of Tank (m³)"]
        y = 5, [description = "Height of the store (m)"]
        h = 1000, [description = "volumetric heat transfer coefﬁcient (W/m³/K)"]
    end

    eqs = [
        Vdot ~ inlet.mdot*ρg
        outlet.mdot ~ inlet.mdot
        Δx ~ y/N
        [x[i] ~ i*Δx for i = 1:N+1]

        Tg[1] ~ inlet.T
        [ρg*Cg*ϵ*D(Tg[i]) ~ -(ρg*Cg*Vdot/A)*((Tg[i+1] - Tg[i-1])/2Δx) + h*(Ts[i] - Tg[i]) - (UA/N)*(Tg[i] - T_amb) for i = 2:N]
        D(Tg[end]) ~ -(ρg*Cg*Vdot/A)*((Tg[end] - Tg[end-1])/Δx) + h*(Ts[end] - Tg[end]) - (UA/N)*(Tg[end] - T_amb)
        [ρs*Cs*(1-ϵ)*D(Ts[i]) ~ h*(Tg[i] - Ts[i]) for i = 1:N+1]
        outlet.T ~ Tg[end]
    ]
    return compose(ODESystem(eqs, t, vars, para;name=name),inlet,outlet)
end

export PackedBed