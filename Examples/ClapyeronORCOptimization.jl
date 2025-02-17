using CarnotCycles, CoolProp, ModelingToolkit, DifferentialEquations, Clapeyron
@independent_variables t
model = cPR(["isopentane","toluene"],idealmodel = ReidIdeal)
load_fluid(model)

start_T = 300;
start_p = 101325
start_mdot = 40 #g/s



@named source = MassSource()
@named compressor = Pump()
@named evap = CarnotCycles.SimpleEvaporator()
@named expander = CarnotCycles.IsentropicExpander()
@named condenser = CarnotCycles.SimpleCondensor()
@named sink = MassSink()

eqs = [
    connect(source.port,compressor.inport)
    connect(compressor.outport,evap.inport)
    connect(evap.outport,expander.inport)
    connect(expander.outport,condenser.inport)
    connect(condenser.outport,sink.port)
]
systems = [source,compressor,evap,expander,condenser,sink]
@named orc = ODESystem(eqs, t, systems=systems)
sys = structural_simplify(orc)
function ORC(x,p)
    @show x
    x_ = x[1]
    z_ = CarnotCycles.mass_to_moles(model,x_,start_mdot)
    start_h = enthalpy(model,start_p,start_T,z_)
    u0 = []
    para = [source.source_pressure=>start_p, source.source_enthalpy => start_h,source.source_mdot => start_mdot, 
            source.source_x => x_,
            compressor.πc => x[3],compressor.η => p[1],
            expander.η => p[2], expander.πc => compressor.πc, 
            evap.ΔT_sh => x[2],
            condenser.ΔT_sc => condenser.T_sat - start_T]
prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)
if sol[evap.T_out] > 400
    return 0.0
end
@show sol[evap.T_out], sol[expander.T_out], sol[condenser.T_out], sol[condenser.T_sat]

@show sol[expander.P]
η_cycle = (sol[expander.P] - sol[compressor.P])/sol[evap.Qdot]
@show η_cycle
return η_cycle
end

x = [0.5,6,37.0]
p = [0.7,0.7]
ORC([0.99,2,8],p)
# using Optimization, OptimizationMetaheuristics
# f = OptimizationFunction(ORC)
# prob = Optimization.OptimizationProblem(f, x, p, lb = [0.6,2.0,1.0], ub = [0.99,30,5.0])
# sol = solve(prob, DE(),use_initial=false, maxiters = 500, maxtime = 1000.0,x_tol = 1e-2,
#         f_tol = 1e-2,f_tol_rel = 1e-2)
