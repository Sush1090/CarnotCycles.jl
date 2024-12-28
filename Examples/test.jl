using CarnotCycles, ModelingToolkit, DifferentialEquations, Clapeyron

@independent_variables t
model = cPR(["isopentane","toluene"],idealmodel = ReidIdeal);
load_fluid(model)

@named src = MassSource()
@named comp = IsentropicCompressor()
@named evap = SimpleEvaporator()
@named sink = MassSink()
eqs = [
    connect(src.port,comp.inport)
    connect(comp.outport,evap.inport)
    connect(evap.outport,sink.port)
]

system=[src,comp,evap,sink]
@named cycle = ODESystem(eqs,t,systems= system)
sys = structural_simplify(cycle)
# para = [sys.src.source_x=>1.0,sys.comp.πc => 2,sys.src.source_mdot => 20]
z_ = CarnotCycles.mass_to_moles(model,[0.6,0.4],100)
h_ = enthalpy(model,101325,300,z_)
para = [sys.src.source_pressure => 101325,sys.src.source_enthalpy => h_, sys.src.source_x => 0.6,
        sys.evap.ΔT_sh => 3, sys.src.source_mdot => 100,sys.comp.πc => 5]
u0 = []

prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)