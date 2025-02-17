
pt_entropy(model::EoSModel,p,T,z) = Clapeyron.entropy(model::EoSModel,p,T,z,phase = "unknown")
@register_symbolic pt_entropy(model::EoSModel,p,T,z)

pt_enthalpy(model::EoSModel,p,T,z) = Clapeyron.enthalpy(model::EoSModel,p,T,z,phase = "unknown")
@register_symbolic pt_enthalpy(model::EoSModel,p,T,z)

pt_mass_density(model::EoSModel,p,T,z) = Clapeyron.mass_density(model::EoSModel,p,T,z,phase = "unknown")
@register_symbolic pt_mass_density(model::EoSModel,p,T,z)

ph_temperature(model::EoSModel,p,h,z) = PH.temperature(model::EoSModel,p,h,z,phase = "unknown")
@register_symbolic ph_temperature(model::EoSModel,p,h,z)

ph_entropy(model::EoSModel,p,h,z) = PH.entropy(model::EoSModel,p,h,z,phase = "unkown")
@register_symbolic ph_entropy(model::EoSModel,p,h,z)

ph_mass_density(model::EoSModel,p,h,z) = PH.mass_density(model::EoSModel,p,h,z,phase = "unkown")
@register_symbolic ph_mass_density(model::EoSModel,p,h,z)

ph_volume(model::EoSModel,p,h,z) = PH.volume(model::EoSModel,p,h,z,phase = "unkown")
@register_symbolic ph_volume(model::EoSModel,p,h,z)

ps_temperature(model::EoSModel,p,s,z) = PS.temperature(model::EoSModel,p,s,z,phase = "unknown")
@register_symbolic ps_temperature(model::EoSModel,p,s,z)

ps_enthalpy(model::EoSModel,p,s,z) = PS.enthalpy(model::EoSModel,p,s,z,phase = "unknown")
@register_symbolic ps_enthalpy(model::EoSModel,p,s,z)

Tproperty_S(model::EoSModel,p,s,z) = Clapeyron.Tproperty(model::EoSModel,p,s,z,entropy,phase = "unkown")
@register_symbolic Tproperty_S(model::EoSModel,p,s,z)

Tproperty_H(model::EoSModel,p,h,z) = Clapeyron.Tproperty(model::EoSModel,p,h,z,enthalpy,phase = "unkown",verbose = false)
@register_symbolic Tproperty_H(model::EoSModel,p,h,z)

Bubble_temperature(model::EoSModel,p,z) = Clapeyron.bubble_temperature(model,p,z,FugBubbleTemperature())[1]
@register_symbolic Bubble_temperature(model::EoSModel,p,z)

Dew_temperature(model::EoSModel,p,z) = Clapeyron.dew_temperature(model,p,z,FugDewTemperature())[1]
@register_symbolic Dew_temperature(model::EoSModel,p,z)

Bubble_pressure(model::EoSModel,T,z) = Clapeyron.bubble_pressure(model,T,z,FugBubblePressure())[1]
@register_symbolic Bubble_pressure(model::EoSModel,T,z)

Dew_pressure(model::EoSModel,T,z) = Clapeyron.dew_pressure(model,T,z,FugDewPressure())[1]
@register_symbolic Dew_pressure(model::EoSModel,T,z)




function CriticalTemperature(model::EoSModel,z)
    if size(z,1) == 1
        return Clapeyron.crit_pure(model)[1]
    end
    if size(z,1) != 1
        return Clapeyron.crit_mix(model,z)[1]
    end
end
@register_symbolic CriticalTemperature(model::EoSModel,z)

function CriticalPressure(model::EoSModel,z)
    if size(z,1) == 1
        return Clapeyron.crit_pure(model)[2]
    end
    if size(z,1) != 1
        return Clapeyron.crit_mix(model,z)[2]
    end
end
@register_symbolic CriticalPressure(model::EoSModel,z)


function LiquidPhaseChecker(model::EoSModel,p,h,z)
    # phase = ClapyeronPhase_PT(model,p,T,z)
    phase = PhaseIdentification(model,p,h,z)
    @assert phase == :liquid "Phase not Liquid - should be liquid here"
    return 1
end
@register_symbolic LiquidPhaseChecker(model::EoSModel,p,h,z)


function GasPhaseChecker(model::EoSModel,p,h,z)
    # phase = ClapyeronPhase_PT(model,p,T,z)
    phase = PhaseIdentification(model,p,h,z)
    @assert phase == :vapour "Phase not gas - should be gas here"
    return 1
end
@register_symbolic GasPhaseChecker(model::EoSModel,p,h,z)

function TwoPhaseChecker(model::EoSModel,p,h,z)
    # phase = ClapyeronPhase_PT(model,p,T,z)
    phase = PhaseIdentification(model,p,h,z)
    @assert phase == :TwoPhase "Phase not Twophase - should be Twophase here"
    return 1
end
@register_symbolic TwoPhaseChecker(model::EoSModel,p,h,z)


function PhaseIdentification(model::EoSModel,p,h,z)
    Tcrit = CriticalTemperature(model,z)
    T = ph_temperature(model,p,h,z)
    
    if (T>=Tcrit)
        return :vapour
    end
    if (T< Tcrit)
        bt = Bubble_temperature(model,p,z)
        dt = Dew_temperature(model,p,z)
        if bt<= T <= dt
            return :TwoPhase
        end
        if dt< T
            return :vapour
        end
        if bt> T
            return :liquid
        end
    end
end
@register_symbolic PhaseIdentification(model::EoSModel,p,h,z)




function qp_enthalpy(model::EoSModel,q,p,z)
    @assert 0<=q "Vapour fraction less than 0. Should be between [0,1]"
    @assert q<=1 "Vapour fraction greater than 1. Should be between [0,1]"
    res = Clapeyron.qp_flash(model,q,p,z)
    h = Clapeyron.enthalpy(model,res)
    return h
end
@register_symbolic qp_enthalpy(model::EoSModel,q,p,z)

function qp_entropy(model::EoSModel,q,p,z)
    @assert 0<=q "Vapour fraction less than 0. Should be between [0,1]"
    @assert q<=1 "Vapour fraction greater than 1. Should be between [0,1]"
    res = Clapeyron.qp_flash(model,q,p,z)
    s = Clapeyron.entropy(model,res)
    return s
end
@register_symbolic qp_entropy(model::EoSModel,q,p,z)


function qp_temperature(model::EoSModel,q,p,z)
    @assert 0<=q "Vapour fraction less than 0. Should be between [0,1]"
    @assert q<=1 "Vapour fraction greater than 1. Should be between [0,1]"
    res = Clapeyron.qp_flash(model,q,p,z)
    T = Clapeyron.temperature(model,res)
    return T
end
@register_symbolic qp_temperature(model::EoSModel,q,p,z)

function qp_mass_density(model::EoSModel,q,p,z)
    @assert 0<=q "Vapour fraction less than 0. Should be between [0,1]"
    @assert q<=1 "Vapour fraction greater than 1. Should be between [0,1]"
    res = Clapeyron.qp_flash(model,q,p,z)
    ρ = Clapeyron.mass_density(model,res)
    return ρ
end
@register_symbolic qp_mass_density(model::EoSModel,q,p,z)

