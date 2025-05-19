abstract type HEXModel end

struct ϵNTU <:HEXModel
    flow::Symbol
    U::Float64
    A::Float64
    function ϵNTU(;flow = :counterflow, U = 0.0, A = 0.0)
        if flow != :counterflow && flow != :parallelflow && flow != :crossflow
            throw(ArgumentError("flow must be :counterflow, :parallelflow, or :crossflow. Check the flow type"))
        end
        new(flow,U,A)
    end
end