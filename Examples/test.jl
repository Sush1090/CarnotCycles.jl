using CarnotCycles, ModelingToolkit, DifferentialEquations, Clapeyron

@independent_variables t
model = cPR(["ethane","methane"],idealmodel = ReidIdeal);
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
para = [sys.src.source_x=>0.5,sys.comp.πc => 4]
u0 = []

prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)