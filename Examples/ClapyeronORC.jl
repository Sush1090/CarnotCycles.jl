using CarnotCycles, CoolProp, ModelingToolkit, DifferentialEquations, Clapeyron
@independent_variables t
model = cPR(["Pentane","toluene"],idealmodel = ReidIdeal)
load_fluid(model)

start_T = 300;
start_p = 101325
start_mdot = 40 #g/s
x_ = 0.9
z_ = CarnotCycles.mass_to_moles(model,[x_,1-x_],start_mdot)
start_h = enthalpy(model,start_p,start_T,z_)


@named source = MassSource()
@named compressor = Pump()
@named evap = CarnotCycles.SimpleEvaporator()
@named expander = CarnotCycles.IsentropicExpander()
@named sink = MassSink()

eqs = [
    connect(source.port,compressor.inport)
    connect(compressor.outport,evap.inport)
    connect(evap.outport,expander.inport)
    connect(expander.outport,sink.port)
]
systems = [source,compressor,evap,expander,sink]
@named orc = ODESystem(eqs, t, systems=systems)
u0 = []
para = [source.source_pressure=>start_p, source.source_enthalpy => start_h,source.source_mdot => start_mdot, 
        source.source_x => x_,
        compressor.πc => 37.0,compressor.η => 0.7,
        expander.η => 0.9, expander.πc => compressor.πc, 
        evap.ΔT_sh => 6]
sys = structural_simplify(orc)
prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)

η_cycle = (sol[expander.P] - sol[compressor.P])/sol[evap.Qdot]

# bool1 = isapprox(sol[compressor.s_in],sol[compressor.s_out])
# bool2 = isapprox(sol[expander.s_in],sol[expander.s_out])