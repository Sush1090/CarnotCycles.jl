
"""
`Compute_cycle_error(p_source,h_source,p_sink,h_sink;reltol = 1e-8)`
    * Computes cycle end point mismatch in states - between sink and source.

    returns `nothing`

    Shows error incase of mismatch
"""
function Compute_cycle_error(p_source,h_source,p_sink,h_sink;reltol = 1e-8)
    err_p = abs(p_source - p_sink)/p_source
    err_h = abs(h_source - h_sink)/h_source
    if err_p > reltol && err_h <reltol
        @warn "Pressure mismatch between cycle source and sink"
    end
    if err_h > reltol && err_p < reltol
        @warn "Enthalpy mismatch between cycle source and sink"
    end
    if err_h > reltol && err_p > reltol
        @warn "Enthalpy and pressure mismatch between cycle source and sink"
    end
    return nothing
end

"""
`Compute_cycle_error(sol::ODESolution,system::Vector{System};reltol = 1e-8)`
"""
function Compute_cycle_error(sol::ODESolution,system::Vector{System};reltol = 1e-8)
    source = system[1]; sink = system[end]
    h_source = sol[source.h][1]; p_source = sol[sink.p][1];
    h_sink = sol[sink.h][1]; p_sink = sol[sink.p][1];
    Compute_cycle_error(p_source,h_source,p_sink,h_sink,reltol=reltol)
end

export Compute_cycle_error


