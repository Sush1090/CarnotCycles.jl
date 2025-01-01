using CarnotCycles, CoolProp, ModelingToolkit, DifferentialEquations, Plots

load_fluid("R134A")
@independent_variables t
T_ = 290; p_ = 101325; N = 100
h_ = PropsSI("H","T",T_,"P",p_,"R134A")
@named source = MassSource()
@named compressor = CarnotCycles.IsentropicCompressor()
@named condensor = CarnotCycles.SimpleCondensor()
@named valve = Valve()
@named evaporator = SimpleEvaporator()
@named sink = MassSink()
@named h2s = Heat2Storage()
@named store = PackedBed(N = N)

systems = [source, compressor,condensor,valve,evaporator,h2s,store,sink]
eqs = [
        connect(source.port,compressor.inport)
        connect(compressor.outport,condensor.inport)
        connect(condensor.outport,valve.inport)
        connect(valve.outport,evaporator.inport)
        connect(evaporator.outport,sink.port)
        connect(condensor.heatport,h2s.heatport)
        connect(h2s.storeport,store.storeport)
]

@named CB_charging = ODESystem(eqs,t,systems = systems)
sys = structural_simplify(CB_charging)
para = [
    source.source_pressure => p_, source.source_enthalpy => h_, source.source_mdot => 0.02, 
    compressor.πc => 5, compressor.η => 0.9,
    condensor.ΔT_sc => 3,
    valve.πc => compressor.πc,
    evaporator.ΔT_sh => 10, 
    h2s.Cp => 4200, h2s.T_in => condensor.T_out - 2, h2s.T_out => condensor.T_in - 2,
    store.Cg => h2s.Cp, store.ρg => 1000
]
Tg_ = ones(N+1).*285; Ts_ = ones(N+1).* 285
u0 = [
    store.Tg => Tg_ , store.Ts => Ts_
]
tspan = (0,5000)
prob = ODEProblem(sys,u0,tspan,para,symbolic_u0 = true)
sol = solve(prob, Rodas4(autodiff=false))

plot(sol[store.x][end],sol[store.Ts][end],label = "solid at charge end",xlabel = "Axis of the storage (m)" , ylabel = "Temperature (K)")
plot!(sol[store.x][end],sol[store.Tg][end],label = "gas at charge end")