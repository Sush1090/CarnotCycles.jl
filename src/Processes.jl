
using CoolProp, ModelingToolkit

"""
`IsentropicCompression(πc, h_in, p_in,fluid,η)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid
    5. `η`    : Isentropic Efficiency

* Output -> Outlet enthalpy after isentropic compression
"""
function IsentropicCompression(πc, h_in, p_in,fluid,η)
    @assert η <= 1 "Efficiency more than 1"
    s_in = PropsSI("S", "H", h_in, "P", p_in, fluid)
    h_is = PropsSI("H", "S", s_in, "P",πc*p_in, fluid)
    h_out = h_in  + (h_is -h_in)/η
    return h_out
end
@register_symbolic IsentropicCompression(πc, h_in, p_in,fluid::String,η)
export IsentropicCompression

function IsentropicCompressionClapeyron(πc, h_in, p_in,z,fluid::EoSModel,η)
    @assert η <= 1 "Efficiency more than 1"
    s_in =  ph_entropy(fluid,p_in,h_in,z)
    h_is =  ps_enthalpy(fluid,p_in*πc,s_in,z)
    h_out = h_in  + (h_is -h_in)/η
    return h_out
end
@register_symbolic IsentropicCompressionClapeyron(πc, h_in, p_in,z,fluid::EoSModel,η)
export IsentropicCompressionClapeyron
"""
`IsentropicExpansion(πc, h_in, p_in,fluid,η)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid
    5. `η`    : Isentropic Efficiency

* Output -> Outlet enthalpy after isentropic expansion
"""
function IsentropicExpansion(πc, h_in, p_in,fluid,η)
    @assert η <= 1 "Efficiency more than 1"
    s_in = PropsSI("S", "H", h_in, "P", p_in, fluid)
    h_is = PropsSI("H", "S", s_in, "P",p_in/πc, fluid)
    h_out = h_in  - η*(h_in - h_is)
    return h_out
end
@register_symbolic IsentropicExpansion(πc, h_in, p_in,fluid::AbstractString,η)
export IsentropicExpansion

function IsentropicExpansionClapeyron(πc, h_in, p_in,z,fluid::EoSModel,η)
    @assert η <= 1 "Efficiency more than 1"
    s_in =  ph_entropy(fluid,p_in,h_in,z)
    h_is =  ps_enthalpy(fluid,p_in/πc,s_in,z)
    h_out = h_in  - η*(h_in - h_is)
    return h_out
end
@register_symbolic IsentropicExpansionClapeyron(πc, h_in, p_in,z,fluid::EoSModel,η)
export IsentropicExpansionClapeyron


"""
`IsochoricCompression(πc, h_in, p_in,fluid)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid

* Output -> Outlet enthalpy after isochoric compression
"""
function IsochoricCompression(πc, h_in, p_in,fluid)
    v_in = 1/PropsSI("D", "H", h_in, "P", p_in, fluid)
    h_out =  PropsSI("H", "D", 1/v_in, "P", πc*p_in, fluid)
    return h_out
end
@register_symbolic IsochoricCompression(πc, h_in, p_in,fluid::AbstractString)
export IsochoricCompression

function IsochoricCompressionClapeyron(πc, h_in, p_in,z::Array,fluid::EoSModel)
    v_in =  ph_volume(fluid,p_in,h_in,z)
    p_out = p_in*πc
    f(h) = ph_volume(fluid,p_out,h,z) - v_in
    prob = Roots.ZeroProblem(f,h_in)
    sol = Roots.solve(prob)
    return sol
end
@register_symbolic IsochoricCompressionClapeyron(πc, h_in, p_in,z::Array,fluid::EoSModel)
export IsochoricCompressionClapeyron

"""
`IsochoricExpansion(πc, h_in, p_in,fluid)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid

* Output -> Outlet enthalpy after isochoric expansion
"""
function IsochoricExpansion(πc, h_in, p_in,fluid)
    v_in = 1/PropsSI("D", "H", h_in, "P", p_in, fluid)
    h_out =  PropsSI("H", "D", 1/v_in, "P", p_in/πc, fluid)
    return h_out
end
@register_symbolic IsochoricExpansion(πc, h_in, p_in,fluid::AbstractString)
export IsochoricExpansion
function IsochoricExpansionClapeyron(πc, h_in, p_in,z::Array,fluid::EoSModel)
    v_in =  ph_volume(fluid,p_in,h_in,z)
    p_out = p_in/πc
    f(h) = ph_volume(fluid,p_out,h,z) - v_in
    prob = Roots.ZeroProblem(f,h_in)
    sol = Roots.solve(prob)
    return sol
end
@register_symbolic IsochoricExpansionClapeyron(πc, h_in, p_in,z::Array,fluid::EoSModel)
export IsochoricExpansionClapeyron

"""
`IsothermalCompression(πc, h_in, p_in,fluid)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid

* Output -> Outlet enthalpy after Isothermal Compression
"""
function IsothermalCompression(πc, h_in, p_in,fluid)
    T_in = PropsSI("T", "H", h_in, "P", p_in, fluid)
    h_out = PropsSI("H", "T", T_in, "P",πc*p_in, fluid)
    return h_out
end
@register_symbolic IsothermalCompression(πc, h_in, p_in,fluid::AbstractString)
export IsothermalCompression

function IsothermalCompressionClapeyron(πc, h_in, p_in,z,fluid::EoSModel)
    T_in =  ph_temperature(fluid,p_in,h_in,z)
    h_out = pt_enthalpy(fluid,p_in*πc,T_in,z)
    return h_out
end
@register_symbolic IsothermalCompressionClapeyron(πc, h_in, p_in,z,fluid::EoSModel)
export IsothermalCompressionClapeyron

function IsothermalExpansionClapeyron(πc, h_in, p_in,z,fluid::EoSModel)
    T_in =  ph_temperature(fluid,p_in,h_in,z)
    h_out = pt_enthalpy(fluid,p_in/πc,T_in,z)
    return h_out
end
@register_symbolic IsothermalExpansionClapeyron(πc, h_in, p_in,z,fluid::EoSModel)
export IsothermalExpansionClapeyron


"""
`IsothermalExpansion(πc, h_in, p_in,fluid)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid

* Output -> Outlet enthalpy after Isothermal Expansion
"""
function IsothermalExpansion(πc, h_in, p_in,fluid)
    T_in = PropsSI("T", "H", h_in, "P", p_in, fluid)
    h_out = PropsSI("H", "T", T_in, "P", p_in/πc, fluid)
    return h_out
end
@register_symbolic IsothermalExpansion(πc, h_in, p_in,fluid::AbstractString)
export IsothermalExpansion


