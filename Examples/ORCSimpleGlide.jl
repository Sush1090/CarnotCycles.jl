using CarnotCycles, ModelingToolkit, Clapeyron, DifferentialEquations


@independent_variables t

fluid = cPR(["isopentane","isobutane"],idealmodel = ReidIdeal)
load_fluid(fluid)

@named source = MassSource()
@named pump = Pump()
@named evaporator = CarnotCycles.SimpleEvaporatorGlide()
@named turbine = IsentropicExpander()
@named condensor = SimpleCondensorGlide()
@named sink = MassSink()

eqs = [
    connect(source.port,pump.inport)
    connect(pump.outport,evaporator.inport)
    connect(evaporator.outport,turbine.inport)
    connect(turbine.outport,condensor.inport)
    connect(condensor.outport,sink.port)
]

systems = [source,pump,evaporator,turbine,condensor,sink]
@named system = ODESystem(eqs,t,systems = systems)
@time sys = structural_simplify(system)


"""
x = [p_source,x_source,πc,ΔT_sh,ΔT_sc]
p = [η_exp,η_pump,T_htf_evap,T_htf_cond,source_mdot]
"""
function ORC(x,p)
try
    T_htf_cond = p[4];
    T_htf_evap = p[3];
    @assert T_htf_evap[1] > T_htf_evap[2]
    @assert T_htf_cond[1] < T_htf_cond[2]
    z_source = CarnotCycles.mass_to_moles(fluid,x[2],p[5])
    T_source = CarnotCycles.Bubble_temperature(fluid,x[1],z_source) - x[5]
    h_source = CarnotCycles.pt_enthalpy(fluid,x[1],T_source,z_source)
    u0 = []
    para = [source.source_pressure => x[1], source.source_enthalpy => h_source, source.source_mdot => p[5], source.source_x => x[2],
        pump.πc => x[3], pump.η => p[2],
        evaporator.ΔT_sh => x[4], evaporator.T_htf_in => T_htf_evap[1], evaporator.T_htf_out => T_htf_evap[2],
        turbine.η => p[1], turbine.πc => pump.πc,
        condensor.T_htf_in => T_htf_cond[1], condensor.T_htf_out => T_htf_cond[2], condensor.ΔT_sc => x[5]
        ]

    prob = SteadyStateProblem(sys,u0,para)
    sol =solve(prob)
    
    # Check inlet of pump to be liquid
    try
        sol[pump.LiquidPhase]
    catch
        return 1e4
    end

    # Check if the temperature profiles in the evaporator and condenser are feasible

    if sol[evaporator.is_feas] == false
        return 1e4
    end

    if sol[condensor.is_feas] == false
        return 1e4
    end
    return @show (sol[turbine.P] + sol[pump.P])/sol[evaporator.Qdot]
catch
    return 1e4
end
end


"""
x = [p_source,x_source,πc,ΔT_sh,ΔT_sc]
p = [η_exp,η_pump,T_htf_evap,T_htf_cond,source_mdot]
"""


x0 = [101325*5,0.5,3,3,3]
para = [0.7,0.7,(390,380),(290,300),20]
lb = [101325,0.01,2,2,2]
ub = [101325*5,0.99,15,10,10]

using Optimization, OptimizationMetaheuristics

fopt = OptimizationFunction(ORC)
prob_opt = OptimizationProblem(fopt,x0,para,lb=lb,ub=ub)
@show sol = solve(prob_opt, ABC(N = 100), maxiters = 100000, maxtime = 900.0)