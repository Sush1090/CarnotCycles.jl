using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
load_fluid("Argon")
start_T =  280; # Temperature at source 
start_mdot = 0.2 #kg/s

source_h = PropsSI("H","P",101325,"T",280,"Argon")

@named source = MassSource()
@named comp = IsentropicCompressor()
@named superheater = IsobaricHeatSource()
@named expander = IsentropicExpander()
@named sink = MassSink()

# Define equations
eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,superheater.inport)
    connect(superheater.outport,expander.inport)
    connect(expander.outport,sink.port)
]
system=[source,comp,superheater,expander,sink] # Define system

@named brayton_cycle = ODESystem(eqs, t, systems=system)
para = [
    source.source_pressure => 101325,source.source_mdot => start_mdot,source.source_enthalpy => source_h,
    comp.πc => 2, comp.η => 0.7,
    expander.πc => comp.πc, expander.η => comp.η, 
   superheater.Qdot=>2e5
]
u0 = []
sys = structural_simplify(brayton_cycle)
prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)


#compute Efficiency of the cycle.
#Note: sign convetion: Power supplied to the system is +ve while from thee system is -ve
@show η = (sol[expander.P].+ sol[comp.P])./sol[superheater.P]

# #Check if the final state is close to the inital state. 
# Compute_cycle_error(sol,system,reltol = 1e-2)
# CarnotCycles.CyclePlot(PhasePlotType_PH(),sol,system)