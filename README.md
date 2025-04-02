# CarnotCycles.jl

[![Build Status](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sush1090.github.io/CarnotCycles.jl/dev/)

This package combines [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) "Acausal Modeling" with [CoolProp.jl](https://github.com/CoolProp/CoolProp.jl) and [Clapyeron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl) to model thermodynamic cycles.


![CarnotCycle_logo](docs/src/Images/CarnotCycles_logo.png)

It is based on acausal framework so users can add their own components for simulation.
```julia
function MyComp(type::abc;name,...)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        MyParas ...
    end
    vars = @variables begin
        MyVars ...
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~ eq1 ...
            outport.h ~ eq2 ...
            ..
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end
```

## Installation

In the Julia prompt, first type `]` and then:

```julia
pkg> add CarnotCycles
```

For the developer version type:

```julia
pkg> add https://github.com/Sush1090/CarnotCycles.jl.git
```



