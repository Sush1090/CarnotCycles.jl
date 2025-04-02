# Cycle Modeling

## Carnot Cycle
As the name of the package is CarnotCycles.jl we would like to show the first as the cycle proposed by Carnot called the [Carnot Cycle](https://en.wikipedia.org/wiki/Carnot_cycle). 

His cycle follows a isothermal exapansion of the gas, isentropic expansion, isothermal compression , and finally isentropic compression.

So we will use Clapeyron.jl for our gas model. Here we choose the gas to be Argon.

```julia
using CarnotCycles, ModelingToolkit, Clapeyron, DifferentialEquations

fluid = cPR(["Argon"],idealmodel = ReidIdeal)
load_fluid(fluid)
@independent_variables t
```

The we choose our processes as components and connect them. A `source` and `sink` is recommended to initiate and close the cycle.
```julia
@named source = MassSource()
@named isothermal_comp =  IsothermalCompressor()
@named isentropic_comp = IsentropicCompressor()
@named isothermal_exp = IsothermalExpander()
@named isentropic_exp = IsentropicExpander()
@named sink = MassSink()

eqs = [
    connect(source.port,isothermal_exp.inport)
    connect(isothermal_exp.outport,isentropic_exp.inport)
    connect(isentropic_exp.outport,isothermal_comp.inport)
    connect(isothermal_comp.outport,isentropic_comp.inport)
    connect(isentropic_comp.outport,sink.port)
]

systems = [source,isothermal_comp,isothermal_exp,isentropic_comp,isentropic_exp,sink]

@named CarnotCycle = ODESystem(eqs, t, systems=systems)
@time sys = structural_simplify(CarnotCycle)
```

Now we state the point at `source`
```julia
πc_1 = 5; πc_2 = 6
source_mdot = 30 #g/s
z_source = CarnotCycles.mass_to_moles(fluid,1,source_mdot)
source_temp = 600; source_pressure = 101325*30; source_h = CarnotCycles.pt_enthalpy(fluid,source_pressure,source_temp,z_source)

para = [
    source.source_enthalpy => source_h, source.source_mdot => source_mdot, source.source_x => 1, source.source_pressure => source_pressure,
    isothermal_exp.πc => πc_1,
    isentropic_exp.πc => πc_2, isentropic_exp.η => 1,
    isothermal_comp.πc => πc_1,
    isentropic_comp.πc => πc_2, isentropic_comp.η=>1
]
u0 = []
prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)
```
## Vapour Compression Cycle

![Simple_VCC](Images/SimpleHP.jpg) 

We first start by loading the adequate packages and the fluid:
```julia
using CarnotCycles, CoolProp, ModelingToolkit, DifferentialEquations
@independent_variables t
load_fluid("R134A")
```


Then we define the source thermodynamic state -  the starting point of the cycle.
```julia
ΔT_sh = 5
p_ = 101325*5; T_ = PropsSI("T","Q",1,"P",p_,"R134A") + ΔT_sh
h_ = PropsSI("H","T",T_,"P",p_,"R134A")
```

The we choose the adequate components for the vapour compression cycle:
```julia
@named source = MassSource()
@named compressor = CarnotCycles.IsentropicCompressor()
@named condensor = CarnotCycles.SimpleCondensor()
@named valve = Valve()
@named evaporator = SimpleEvaporator()
@named sink = MassSink()
```

Then we connect them in necessary order:

```julia
systems = [source, compressor,condensor,valve,evaporator,sink]
eqs = [
        connect(source.port,compressor.inport)
        connect(compressor.outport,condensor.inport)
        connect(condensor.outport,valve.inport)
        connect(valve.outport,evaporator.inport)
        connect(evaporator.outport,sink.port)
]

@named VCC = ODESystem(eqs,t,systems = systems)
sys = structural_simplify(VCC)
```

Then we choose the parameters of the system:

```julia
para = [
    source.source_pressure => p_, source.source_enthalpy => h_, source.source_mdot => 0.02, 
    compressor.πc => 3, compressor.η => 0.7,
    condensor.ΔT_sc => 3,
    valve.πc => compressor.πc,
    evaporator.ΔT_sh => ΔT_sh, 
]
```

Then we proceed to solve the problem:

```julia
u0 = []
prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)
```

To get the Coeffecient of Performace of the cycle: 
```julia
julia> COP = sol[condensor.Qdot]/sol[compressor.P]
-5.096928812859646
```

---
**NOTE**

Energy given to the fluid is +ve while given by the fluid is -ve. Hence the COP is negative

---


### Plotting the Cycle


To model a cycle with Claperyon.jl only change the fluid
Example: 
```julia
using CarnotCycles, ModelingToolkit, DifferentialEquations, Clapeyron
@independent_variables t
model = cPR(["Pentane","toluene"],idealmodel = ReidIdeal)
load_fluid(model)
```

