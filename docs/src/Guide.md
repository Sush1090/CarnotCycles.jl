# Guide

## Loading Fluids
Two type of fluid model backends are supported - Clapeyron.jl and CoolProp.jl.

To load a [Clapyeron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl) backend fluid model do the following: 
```julia
using CarnotCycles, Clapeyron
fluid = cPR(["isopentane","isobutane"],idealmodel = ReidIdeal)
load_fluid(fluid)
```
As of now we support Clapeyon mixtures of up to 2 components.

For [CoolProp.jl](https://github.com/CoolProp/CoolProp.jl) backend fluid properties simply pass the fluid name as a string as follows:
```julia
using CarnotCycles, CoolProp
load_fluid("R601")
```
---
**NOTE**

Once the fluid model is chosen through the simulation the underlying components are chosen based on the fluid model.

---



## Loading Components
A cycle consists of various components through which the fluid passes. These components are modeled using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).

To load any component follow ModelingToolkit.jl's `@named` framework. For example: 

```julia
using ModelingToolkit, CarnotCycles

@named evaporator = SimpleEvaporator()
```

For detailed understanding of the acausal framework follow the documentation of [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).


## Port Variables
Every component has atleast two fluid ports `inport` and `outport` (except `MassSource` and `MassSink` which has only one, `port`). 

For a fluid port based over CoolProp backend we have the ports variables to be `h`, `p` and `mdot` (enthalpy, pressure, and mass flow rate). 

For a fluid port based over Clapeyron backend we have the ports variables to be `h`, `p`, `mdot`, and `x` (enthalpy, pressure,mass flow rate, and mass fraction of first fluid).

---
**NOTE**

The units for Coolprop fluid for mass flow rate is kg/s while for Clapeyron is g/s. 

Clapeyron works with total values i.e. total enthalpy (J), entropy (J/K) etc.. while CoolProp works with specific values i.e. specific enthalpy (J/kg), specific entropy (J/kg/K)

---

## Plotting

## Making new components

It is based on acausal framework so users can add their own components for simulation.
```julia
function MyComp(;name,...)
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
