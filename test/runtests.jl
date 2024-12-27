using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp, Clapeyron
using Test


@testset "Global fluid setting" begin
    load_fluid("R601")
    @test CarnotCycles.set_fluid isa AbstractString

    model = cPR(["ethane"]) 
    load_fluid(model)
    @test CarnotCycles.set_fluid isa EoSModel
    @test CarnotCycles.Nc == size(model.components,1)

    model2 = cPR(["ethane","methane"])
    load_fluid(model2)
    @test CarnotCycles.set_fluid isa EoSModel
    @test CarnotCycles.Nc == size(model2.components,1)
end

@testset "Processes - CoolProp" begin
    fluid = "R134A"
    T1 = 300; p1 = 101325; 
    s1 =  CoolProp.PropsSI("S","T",T1,"P",p1,fluid)
    h1 =  CoolProp.PropsSI("H","T",T1,"P",p1,fluid)
    πc = 2; 
    h2 = IsentropicCompression(πc,h1,p1,fluid,1.0)
    s2 = CoolProp.PropsSI("S","H",h2,"P",p1*πc,fluid)
    @test isapprox(s1,s2,atol = 1e-5)
    h3 = IsentropicExpansion(πc,h2,p1*πc,fluid,1.0)
    s3 = CoolProp.PropsSI("S","H",h3,"P",p1,fluid)
    @test isapprox(s1,s3,atol = 1e-5)
    h_isothermal = IsothermalCompression(πc,h1,p1,fluid)
    T_isothermal = CoolProp.PropsSI("T","H",h_isothermal,"P",p1*πc,fluid)
    @test isapprox(T_isothermal,T1,atol = 1e-5)

end


@testset "Processes - Clapeyron - Single Component" begin
    fluid = cPR(["ethane"],idealmodel = ReidIdeal)
    z1 = [5.0]
    T1 = 300; p1 = 101325; h1 = enthalpy(fluid,p1,T1,z1);
    πc = 2; s1 = entropy(fluid,p1,T1,z1)
    h2= IsentropicCompressionClapeyron(πc,h1,p1,z1,fluid,1.0)
    s2 = CarnotCycles.ph_entropy(fluid,p1*πc,h2,z1)
    @test isapprox(s1,s2,atol=1e-5)
end


@testset "Processes - Clapeyron - Dual Component" begin
    fluid = cPR(["ethane","methane"],idealmodel = ReidIdeal)
    z1 = [5.0,5.0]
    T1 = 300; p1 = 101325; h1 = enthalpy(fluid,p1,T1,z1);
    πc = 2; s1 = entropy(fluid,p1,T1,z1)
    h2= IsentropicCompressionClapeyron(πc,h1,p1,z1,fluid,1.0)
    s2 = CarnotCycles.ph_entropy(fluid,p1*πc,h2,z1)
    @test isapprox(s1,s2,atol=1e-5)

    h_isothermal = IsothermalCompressionClapeyron(πc,h1,p1,z1,fluid)
    T_isothermal = CarnotCycles.ph_temperature(fluid,p1*πc,h_isothermal,z1)
    @test isapprox(T1,T_isothermal,atol = 1e-5)
end
# @testset "Isentropic Process - CoolProp" begin
#     fluid = "R134A"
#     load_fluid(fluid)
#     @independent_variables t
#     start_T = 300;
#     start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3
#     ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T;
#     @assert ΔT_subcool > 1e-3
#     start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s


#     @named source = MassSource()
#     @named comp = IsentropicCompressor()
#     @named sink = MassSink()

#     eqs = [
#         connect(source.port,comp.inport)
#         connect(comp.outport,exp.inport)
#         connect(exp.outport,sink.port)
#     ]
#     systems = [source,comp,exp,sink]
#     @named test_isentropic = ODESystem(eqs, t, systems=systems)
#     u0 = []
#     para = [source.source_press=>start_p, source.source_enthalpy = start_h,source.source_mdot => start_mdot,comp.]
#     tspan = (0.0, 1.0)
#     sys = structural_simplify(test_isentropic)
#     prob = ODEProblem(sys,u0,tspan)
#     sol = solve(prob)
#     bool1 = isapprox(sol[comp.s_in][1],sol[comp.s_out][1])
#     bool2 = isapprox(sol[exp.s_in][1],sol[exp.s_out][1])
#     @test bool1 == true
#     @test bool2 == true
#     @test isapprox(sol[source.p][1],sol[sink.p][1])
#     @test isapprox(sol[source.h][1],sol[sink.h][1])
#     if _system.η == 1
#         @test isapprox(sol[comp.s_in][1],sol[comp.s_out][1])
#         @test isapprox(sol[exp.s_in][1],sol[exp.s_out][1])
#     end
#     @test isapprox(sol[comp.p_out][1]/sol[comp.p_in][1],_system.πc)
#     @test isapprox(sol[exp.p_in][1]/sol[exp.p_out][1],_system.πc)
# end


# @testset "Evaporator" begin
#     out_phase = "gas"
#     in_phase_liquid = "liquid"
#     ΔT_sh = 5
#     fluid = "R134A"
#     @independent_variables t
#     start_T = 300;
#     start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3
#     start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s

#     @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
#     @named evporator = SimpleEvaporator(ΔT_sh = ΔT_sh,fluid = fluid)
#     @named sink = MassSink(fluid = fluid)

#     eqs = [
#         connect(source.port,evporator.inport)
#         connect(evporator.outport,sink.port)
#     ]
#     systems = [source,evporator,sink]
#     @named test_evap = ODESystem(eqs, t, systems=systems)
#     u0 = []
#     tspan = (0.0, 1.0)
#     sys = structural_simplify(test_evap)
#     prob = ODEProblem(sys,u0,tspan)
#     sol = solve(prob)

#     out_phase_test = PhaseSI("T",sol[evporator.T_out][1],"P",sol[evporator.p_out][1],fluid)
#     in_phase_test = PhaseSI("T",sol[evporator.T_in][1],"P",sol[evporator.p_in][1],fluid)

#     @test out_phase_test == out_phase
#     @test isapprox(sol[evporator.T_out][1]-ΔT_sh,sol[evporator.T_sat][1])
#     @test in_phase_test == in_phase_liquid
#     @test sol[evporator.h_out][1] >= sol[evporator.h_in][1]
# end


# @testset "Condensor" begin
#     out_phase = "liquid"
#     in_phase = "gas"
#     ΔT_sc = 3
#     fluid = "R134A"
#     @independent_variables t
#     start_T = 300;
#     start_p = PropsSI("P","Q",0,"T",start_T,fluid) - 101325
#     start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s

#     @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
#     @named condensor = SimpleCondensor(ΔT_sc = ΔT_sc,fluid = fluid)
#     @named sink = MassSink(fluid = fluid)

#     eqs = [
#         connect(source.port,condensor.inport)
#         connect(condensor.outport,sink.port)
#     ]
#     systems = [source,condensor,sink]
#     @named test_condensor = ODESystem(eqs, t, systems=systems)
#     u0 = []
#     tspan = (0.0, 1.0)
#     sys = structural_simplify(test_condensor)
#     prob = ODEProblem(sys,u0,tspan)
#     sol = solve(prob)

#     out_phase_test = PhaseSI("T",sol[condensor.T_out][1],"P",sol[condensor.p_out][1],fluid)
#     in_phase_test = PhaseSI("T",sol[condensor.T_in][1],"P",sol[condensor.p_in][1],fluid)

#     @test out_phase_test == out_phase
#     @test isapprox(sol[condensor.T_out][1]+ΔT_sc,sol[condensor.T_sat][1])
#     @test in_phase_test == in_phase
#     @test sol[condensor.h_out][1] <= sol[condensor.h_in][1]
# end


# @testset "Valve - Isenthalpic" begin
#     fluid = "R134A"
#     @independent_variables t
#     start_T = 300;
#     start_p = PropsSI("P","Q",0,"T",start_T,fluid) 
#     start_h = PropsSI("H","Q",0,"P",start_p,fluid); start_mdot = 0.2 #kg/s

#     valve_system  = IsenthalpicExpansionValve(4.5)

#     @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
#     @named exp = Valve(valve_system,fluid= fluid)
#     @named sink = MassSink(fluid = fluid)

#     eqs = [
#         connect(source.port,exp.inport)
#         connect(exp.outport,sink.port)
#     ]
#     systems = [source,exp,sink]
#     @named test_condensor = ODESystem(eqs, t, systems=systems)
#     u0 = []
#     tspan = (0.0, 1.0)
#     sys = structural_simplify(test_condensor)
#     prob = ODEProblem(sys,u0,tspan)
#     sol = solve(prob)
#     @test isapprox(sol[source.h][1],sol[sink.h][1]) 
# end

# @testset "Processes-expansion" begin
#     fluid = "R134A"
#     πc = 4
#     p_in = 101325*πc; T_in = 450;
#     h_in = PropsSI("H","T",T_in,"P",p_in,fluid)
#     s_in =  PropsSI("S","T",T_in,"P",p_in,fluid)
#     v_in = 1/PropsSI("D","T",T_in,"P",p_in,fluid)

#     h_out_isen = IsentropicExpansion(πc,h_in,p_in,fluid,1.0)
#     s_out_isen = PropsSI("S","H",h_out_isen,"P",p_in/πc,fluid)
#     @test isapprox(s_out_isen,s_in)
#     h_out_isothermal = IsothermalExpansion(πc,h_in,p_in,fluid)
#     T_out_isothermal = PropsSI("T","H",h_out_isothermal,"P",p_in/πc,fluid)
#     @test isapprox(T_out_isothermal,T_in)
#     h_out_isochoric = IsochoricExpansion(πc,h_in,p_in,fluid)
#     v_out_isothermal = 1/PropsSI("D","H",h_out_isochoric,"P",p_in/πc,fluid)
#     @test isapprox(v_out_isothermal,v_in)

#     fluid = "Argon"
#     πc = 4
#     p_in = 101325*πc; T_in = 450;
#     h_in = PropsSI("H","T",T_in,"P",p_in,fluid)
#     s_in =  PropsSI("S","T",T_in,"P",p_in,fluid)
#     v_in = 1/PropsSI("D","T",T_in,"P",p_in,fluid)

#     h_out_isen = IsentropicExpansion(πc,h_in,p_in,fluid,1.0)
#     s_out_isen = PropsSI("S","H",h_out_isen,"P",p_in/πc,fluid)
#     @test isapprox(s_out_isen,s_in)
#     h_out_isothermal = IsothermalExpansion(πc,h_in,p_in,fluid)
#     T_out_isothermal = PropsSI("T","H",h_out_isothermal,"P",p_in/πc,fluid)
#     @test isapprox(T_out_isothermal,T_in)
#     h_out_isochoric = IsochoricExpansion(πc,h_in,p_in,fluid)
#     v_out_isothermal = 1/PropsSI("D","H",h_out_isochoric,"P",p_in/πc,fluid)
#     @test isapprox(v_out_isothermal,v_in)
# end

# @testset "Processes-compression" begin
#     fluid = "Argon"
#     πc = 4
#     p_in = 101325; T_in = 450;
#     h_in = PropsSI("H","T",T_in,"P",p_in,fluid)
#     s_in =  PropsSI("S","T",T_in,"P",p_in,fluid)
#     v_in = 1/PropsSI("D","T",T_in,"P",p_in,fluid)

#     h_out_isen = IsentropicCompression(πc,h_in,p_in,fluid,1.0)
#     s_out_isen = PropsSI("S","H",h_out_isen,"P",p_in*πc,fluid)
#     @test isapprox(s_out_isen,s_in)
#     h_out_isothermal = IsothermalCompression(πc,h_in,p_in,fluid)
#     T_out_isothermal = PropsSI("T","H",h_out_isothermal,"P",p_in*πc,fluid)
#     @test isapprox(T_out_isothermal,T_in)
#     h_out_isochoric = IsochoricCompression(πc,h_in,p_in,fluid)
#     v_out_isothermal = 1/PropsSI("D","H",h_out_isochoric,"P",p_in*πc,fluid)
#     @test isapprox(v_out_isothermal,v_in)


#     fluid = "R134A"
#     πc = 2
#     p_in = 101325; T_in = 250;
#     h_in = PropsSI("H","T",T_in,"P",p_in,fluid)
#     s_in =  PropsSI("S","T",T_in,"P",p_in,fluid)
#     v_in = 1/PropsSI("D","T",T_in,"P",p_in,fluid)

#     h_out_isen = IsentropicCompression(πc,h_in,p_in,fluid,1.0)
#     s_out_isen = PropsSI("S","H",h_out_isen,"P",p_in*πc,fluid)
#     @test isapprox(s_out_isen,s_in)
#     h_out_isothermal = IsothermalCompression(πc,h_in,p_in,fluid)
#     T_out_isothermal = PropsSI("T","H",h_out_isothermal,"P",p_in*πc,fluid)
#     @test isapprox(T_out_isothermal,T_in)
#     h_out_isochoric = IsochoricCompression(πc,h_in,p_in,fluid)
#     v_out_isothermal = 1/PropsSI("D","H",h_out_isochoric,"P",p_in*πc,fluid)
#     @test isapprox(v_out_isothermal,v_in)
# end

# @testset "Three-faced valve" begin
#     fluid = "R134A"
#     @load_fluid "R134A"
#     @independent_variables t
#     start_T = 300;
#     start_p = PropsSI("P","Q",0,"T",start_T,fluid) 
#     start_h = PropsSI("H","Q",0,"P",start_p,fluid); start_mdot = 0.2 #kg/s

#     @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot)
#     @named three_split = Valve(ThreeFacedValveSplit(ratio = [0.3,0.7]))
#     @named three_combine = Valve(ThreeFacedValveCombine())
#     @named sink = MassSink()

#     eqs = [
#         connect(source.port,three_split.inport)
#         connect(three_split.outport1,three_combine.inport1)
#         connect(three_split.outport2,three_combine.inport2)
#         connect(three_combine.outport,sink.port)
#     ]
#     systems = [source,three_split,three_combine,sink]

#     @named three = ODESystem(eqs, t, systems=systems)
#     u0 = []
#     tspan = (0.0, 1.0)
#     sys = structural_simplify(three)
#     prob = ODEProblem(sys,u0,tspan)
#     sol = solve(prob)

#     @test isapprox(sol[source.p][1],sol[sink.p][1])
#     @test isapprox(sol[source.h][1],sol[sink.h][1])
#     @test isapprox(sol[source.mdot][1],sol[sink.mdot][1])
# end


# @testset "Process component" begin
#     fluid = "Argon"
#     @load_fluid "Argon"
#     _system = Isentropic_η(η = 1,πc = 5)
#     _isochoric = Isochoric_comp(πc = 5)
#     _isothermal = Isothermal_comp(πc =5)
#     @independent_variables t
#     start_T = 300;
#     start_p = 101325
#     start_h = PropsSI("H","T",start_T,"P",start_p,fluid);start_mdot = 0.2

#     @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot)
#     @named comp_isen = Compressor(_system)
#     @named exp_isen = Expander(_system)
#     @named comp_ther = Compressor(_isothermal)
#     @named exp_ther = Expander(_isothermal)
#     @named comp_chor = Compressor(_isochoric)
#     @named exp_chor = Expander(_isochoric)
#     @named sink = MassSink()

#     eqs = [
#         connect(source.port,comp_isen.inport)
#         connect(comp_isen.outport,exp_isen.inport)
#         connect(exp_isen.outport,comp_ther.inport)
#         connect(comp_ther.outport,exp_ther.inport)
#         connect(exp_ther.outport,comp_chor.inport)
#         connect(comp_chor.outport,exp_chor.inport)
#         connect(exp_chor.outport,sink.port)
#     ]
#     systems = [source,comp_isen,exp_isen,comp_ther,exp_ther,comp_chor,exp_chor,sink]
#     @named test_processes = ODESystem(eqs, t, systems=systems)
#     u0 = []
#     tspan = (0.0, 1.0)
#     sys = structural_simplify(test_processes)
#     prob = ODEProblem(sys,u0,tspan)
#     sol = solve(prob)
#     @test isapprox(sol[comp_isen.s_in][1],sol[comp_isen.s_out][1])
#     @test isapprox(sol[exp_isen.s_in][1],sol[exp_isen.s_out][1])
#     @test isapprox(sol[comp_ther.T_in][1],sol[comp_ther.T_out][1])
#     @test isapprox(sol[exp_ther.T_in][1],sol[exp_ther.T_out][1])
#     @test isapprox(sol[comp_chor.ρ_in][1],sol[comp_chor.ρ_out][1])
#     @test isapprox(sol[exp_chor.ρ_in][1],sol[exp_chor.ρ_out][1])
#     @test isapprox(sol[source.h][1],sol[sink.h][1])
#     @test isapprox(sol[source.p][1],sol[sink.p][1])

# end

@testset "MassSource initilizations" begin
    
end