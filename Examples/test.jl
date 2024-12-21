using CarnotCycles, ModelingToolkit, DifferentialEquations, Clapeyron

@independent_variables t
model = cPR(["ethane"],idealmodel = ReidIdeal);
load_fluid(model)

@named src = MassSource()
@named comp = IsentropicCompressor()
@named sink = MassSink()
eqs = [
    connect(src.port,comp.inport)
    connect(comp.outport,sink.port)
]

system=[src,comp,sink]
@named cycle = ODESystem(eqs,t,systems= system)
sys = structural_simplify(cycle)
para = [sys.src.source_x=>0.5,sys.comp.πc => 2,sys.src.source_mdot => 20]
# para = [sys.src.source_pressure => 101325,sys.src.source_enthalpy => 1e5,
#         sys.comp.πc => 3,sys.comp.η=>0.8, sys.src.source_mdot => 0.02]
u0 = []

prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)