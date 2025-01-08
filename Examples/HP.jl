

using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
load_fluid("propane")

T_ = 290; p_ = 101325;
h_ = PropsSI("H","T",T_,"P",p_,"propane")
@named source = MassSource()
@named compressor = CarnotCycles.IsentropicCompressor()
@named condensor = CarnotCycles.SimpleCondensor()
@named valve = Valve()
@named evaporator = SimpleEvaporator()
@named sink = MassSink()

systems = [source, compressor,condensor,valve,evaporator,sink]
eqs = [
        connect(source.port,compressor.inport)
        connect(compressor.outport,condensor.inport)
        connect(condensor.outport,valve.inport)
        connect(valve.outport,evaporator.inport)
        connect(evaporator.outport,sink.port)
]
@named HP = ODESystem(eqs,t,systems = systems)
sys = structural_simplify(HP)

para = [evaporator.ΔT_sh => 5, condensor.ΔT_sc => 5, compressor.πc => 5, valve.πc => 5, source.source_mdot => 0.02, source.source_pressure => p_, source.source_enthalpy => h_]
u0 = []
prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)
COP = sol[condensor.Qdot]/sol[compressor.P]
print("COP : $COP")