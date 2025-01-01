

using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
load_fluid("propane")

T_ = 290; p_ = 101325; N = 2
h_ = PropsSI("H","T",T_,"P",p_,"R134A")
@named source = MassSource()
@named compressor = CarnotCycles.IsentropicCompressor()
@named condensor = CarnotCycles.SimpleCondensor()
@named valve = Valve()
# @named evaporator = SimpleEvaporator()
@named sink = MassSink()

systems = [source, compressor,condensor,valve,sink]#,condensor,valve,evaporator,sink]
eqs = [
        connect(source.port,compressor.inport)
        connect(compressor.outport,condensor.inport)
        connect(condensor.outport,valve.inport)
        connect(valve.outport,evaporator.inport)
        connect(evaporator.outport,sink.port)
]
@named HP = ODESystem(eqs,t,systems = systems)
sys = structural_simplify(HP)