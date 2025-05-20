using CarnotCycles, ModelingToolkit, SteadyStateDiffEq, CoolProp, Clapeyron
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
    h_isothermal = IsothermalCompressionClapeyron(πc,h1,p1,z1,fluid)
    T_isothermal = CarnotCycles.ph_temperature(fluid,p1*πc,h_isothermal,z1)
    @test isapprox(T1,T_isothermal,atol = 1e-5)
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


@testset "Joining Components - CoolProp" begin
    load_fluid("R134A")
    fluid  =  "R134A"
    @independent_variables t
    start_T = 300;
    start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3
    ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T;
    @assert ΔT_subcool > 1e-3
    # start_h = PropsSI("H","T",start_T,"P",start_p,fluid); 
    start_mdot = 0.2 #kg/s


    @named source = MassSource()
    @named compressor = IsentropicCompressor()
    @named expander = CarnotCycles.IsentropicExpander()
    @named sink = MassSink()

    eqs = [
        connect(source.port,compressor.inport)
        connect(compressor.outport,expander.inport)
        connect(expander.outport,sink.port)
    ]
    systems = [source,compressor,expander,sink]
    @named test_isentropic = ODESystem(eqs, t, systems=systems)
    u0 = []
    para = [source.source_pressure=>start_p, source.source_temperature => start_T,source.source_mdot => start_mdot,compressor.πc => 5.0,compressor.η => 1.0,
            expander.η => 1.0, expander.πc => compressor.πc]
    sys = structural_simplify(test_isentropic)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)
    bool1 = isapprox(sol[compressor.s_in],sol[compressor.s_out])
    bool2 = isapprox(sol[expander.s_in],sol[expander.s_out])
    @test bool1 == true
    @test bool2 == true
    @test isapprox(sol[source.p],sol[sink.p])
    @test isapprox(sol[source.h],sol[sink.h])

    @test isapprox(sol[compressor.s_in],sol[compressor.s_out])
     @test isapprox(sol[expander.s_in],sol[expander.s_out])

end


@testset "Joining Components - Clapyeron - single fluid component- compressor-expander" begin
    model = cPR(["ethane"],idealmodel = ReidIdeal)
    load_fluid(model)
    @independent_variables t
    start_T = 300;
    start_p = 101325
    start_mdot = 20 #g/s
    z_ = CarnotCycles.mass_to_moles(model,1,start_mdot)
    start_h = enthalpy(model,start_p,start_T,z_)


    @named source = MassSource()
    @named compressor = IsentropicCompressor()
    @named expander = CarnotCycles.IsentropicExpander()
    @named sink = MassSink()

    eqs = [
        connect(source.port,compressor.inport)
        connect(compressor.outport,expander.inport)
        connect(expander.outport,sink.port)
    ]
    systems = [source,compressor,expander,sink]
    @named test_isentropic = ODESystem(eqs, t, systems=systems)
    u0 = []
    para = [source.source_pressure=>start_p, source.source_temperature => start_T,source.source_mdot => start_mdot, source.source_x => 1,
            compressor.πc => 5.0, compressor.η => 1.0,
            expander.η => 1.0, expander.πc => compressor.πc]
    sys = structural_simplify(test_isentropic)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)
    bool1 = isapprox(sol[compressor.s_in],sol[compressor.s_out],atol=1e-4)
    bool2 = isapprox(sol[expander.s_in],sol[expander.s_out],atol= 1e-4)
    @test bool1 == true
    @test bool2 == true
    @test isapprox(sol[source.p],sol[sink.p],atol = 1e-4)
    @test isapprox(sol[source.h],sol[sink.h],atol = 1e-4)
end


@testset "Joining Components - Clapyeron - two fluid component -comp-exapander" begin
    model = cPR(["ethane","methane"],idealmodel = ReidIdeal)
    load_fluid(model)
    @independent_variables t
    start_T = 300;
    start_p = 101325
    start_mdot = 20 #g/s
    z_ = CarnotCycles.mass_to_moles(model,0.6,start_mdot)
    start_h = enthalpy(model,start_p,start_T,z_)


    @named source = MassSource()
    @named compressor = IsentropicCompressor()
    @named expander = CarnotCycles.IsentropicExpander()
    @named sink = MassSink()

    eqs = [
        connect(source.port,compressor.inport)
        connect(compressor.outport,expander.inport)
        connect(expander.outport,sink.port)
    ]
    systems = [source,compressor,expander,sink]
    @named test_isentropic = ODESystem(eqs, t, systems=systems)
    u0 = []
    para = [source.source_pressure=>start_p, source.source_temperature => start_T,source.source_mdot => start_mdot,compressor.πc => 5.0,
        compressor.η => 1.0,source.source_x => 0.6,
        expander.η => 1.0, expander.πc => compressor.πc]
    sys = structural_simplify(test_isentropic)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)
    bool1 = isapprox(sol[compressor.s_in],sol[compressor.s_out],atol = 1e-4)
    bool2 = isapprox(sol[expander.s_in],sol[expander.s_out],atol = 1e-4)
    @test bool1 == true
    @test bool2 == true
    @test isapprox(sol[source.p],sol[sink.p],atol = 1e-4)
    @test isapprox(sol[source.h],sol[sink.h],atol = 1e-4)
end



@testset "Glide Evaporator CoolProp" begin
    @independent_variables t
    load_fluid("Water")
    T_start = 370; p = 101325;
    h_start = PropsSI("H","T",T_start,"P",p,"Water")
    start_mdot = 1e-2 #kg/s;
    
    @named source = MassSource()
    @named evap = SimpleEvaporatorGlide()
    @named sink = MassSink()

    eqs = [
        connect(source.port,evap.inport)
        connect(evap.outport,sink.port)
    ]
    systems = [source,evap,sink]
    para = [source.source_pressure => p, source.source_temperature => T_start,source.source_mdot => start_mdot,
    evap.ΔT_sh => 2, evap.T_htf_in => 400, evap.T_htf_out => 380]
    u0 = []
    @named model = ODESystem(eqs, t, systems=systems)
    sys = structural_simplify(model)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    @test sol[evap.is_feas] == true
end


@testset "Glide Evaporator Clapyeron" begin
    @independent_variables t
    fluid = cPR(["Water"],idealmodel = ReidIdeal)
    load_fluid(fluid)
    T = 370; p = 101325;
    start_mdot = 20 #g/s
    z_ = CarnotCycles.mass_to_moles(fluid,1,start_mdot)
    start_h = enthalpy(fluid,p,T,z_)
    
    @named source = MassSource()
    @named evap = SimpleEvaporatorGlide()
    @named sink = MassSink()

    eqs = [
        connect(source.port,evap.inport)
        connect(evap.outport,sink.port)
    ]
    systems = [source,evap,sink]
    para = [source.source_pressure => p, source.source_temperature => T,source.source_mdot => start_mdot, source.source_x => 1,
    evap.ΔT_sh => 2, evap.T_htf_in => 400, evap.T_htf_out => 380]
    u0 = []
    @named model = ODESystem(eqs, t, systems=systems)
    sys = structural_simplify(model)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    @test sol[evap.is_feas] == true
end

@testset "Glide Evaporator Clapeyron two components" begin
    @independent_variables t
    fluid = cPR(["isobutane","isopentane"],idealmodel=ReidIdeal)
    load_fluid(fluid)
    T = 320; p = 101325*5.5 ;
    start_mdot = 20 #g/s
    x_ = 0.7
    z_ = CarnotCycles.mass_to_moles(fluid,x_,start_mdot)
    start_h = enthalpy(fluid,p,T,z_)
    
    @named source = MassSource()
    @named evap = SimpleEvaporatorGlide()
    @named sink = MassSink()

    eqs = [
        connect(source.port,evap.inport)
        connect(evap.outport,sink.port)
    ]
    systems = [source,evap,sink]
    para = [source.source_pressure => p, source.source_temperature => T,source.source_mdot => start_mdot, source.source_x => x_,
    evap.ΔT_sh => 2, evap.T_htf_in => 340, evap.T_htf_out => 330]
    u0 = []
    @named model = ODESystem(eqs, t, systems=systems)
    sys = structural_simplify(model)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    @test sol[evap.is_feas] == true
end


@testset "Glide Condensor CoolProp" begin
    @independent_variables t
    load_fluid("Water")
    T_start = 375; p = 101325;
    h_start = PropsSI("H","T",T_start,"P",p,"Water")
    start_mdot = 1e-2 #kg/s;
    
    @named source = MassSource()
    @named cond = SimpleCondensorGlide()
    @named sink = MassSink()

    eqs = [
        connect(source.port,cond.inport)
        connect(cond.outport,sink.port)
    ]
    systems = [source,cond,sink]
    para = [source.source_pressure => p, source.source_temperature => T_start,source.source_mdot => start_mdot,
    cond.ΔT_sc => 2, cond.T_htf_in => 370, cond.T_htf_out => 372]
    u0 = []
    @named model = ODESystem(eqs, t, systems=systems)
    sys = structural_simplify(model)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    @test sol[cond.is_feas] == true
end


@testset "Glide Condensor Clapyeron" begin
    @independent_variables t
    fluid = cPR(["Water"],idealmodel = ReidIdeal)
    load_fluid(fluid)
    T = 375; p = 101325;
    start_mdot = 20 #g/s
    z_ = CarnotCycles.mass_to_moles(fluid,1,start_mdot)
    start_h = enthalpy(fluid,p,T,z_)
    
    @named source = MassSource()
    @named cond = SimpleCondensorGlide()
    @named sink = MassSink()

    eqs = [
        connect(source.port,cond.inport)
        connect(cond.outport,sink.port)
    ]
    systems = [source,cond,sink]
    para = [source.source_pressure => p, source.source_temperature => T,source.source_mdot => start_mdot, source.source_x => 1,
    cond.ΔT_sc => 2, cond.T_htf_in => 370, cond.T_htf_out => 372]
    u0 = []
    @named model = ODESystem(eqs, t, systems=systems)
    sys = structural_simplify(model)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    @test sol[cond.is_feas] == true
end

@testset "Glide Condensor Clapeyron two components" begin
    @independent_variables t
    fluid = cPR(["isobutane","isopentane"],idealmodel=ReidIdeal)
    load_fluid(fluid)
    T = 335; p = 101325*5.5 ;
    start_mdot = 20 #g/s
    x_ = 0.7
    z_ = CarnotCycles.mass_to_moles(fluid,x_,start_mdot)
    start_h = enthalpy(fluid,p,T,z_)
    
    @named source = MassSource()
    @named cond = CarnotCycles.SimpleCondensorGlide()
    @named sink = MassSink()

    eqs = [
        connect(source.port,cond.inport)
        connect(cond.outport,sink.port)
    ]
    systems = [source,cond,sink]
    para = [source.source_pressure => p, source.source_temperature => T,source.source_mdot => start_mdot, source.source_x => x_,
    cond.ΔT_sc => 2, cond.T_htf_in => 315, cond.T_htf_out => 330]
    u0 = []
    @named model = ODESystem(eqs, t, systems=systems)
    sys = structural_simplify(model)
    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    @test sol[cond.is_feas] == true
end

@testset "Pipes" begin


end

@testset "eNTU" begin

@test CarnotCycles.eNTU(300,0.9,:counterflow) ≈ 1.0 atol = 1e-5
@test CarnotCycles.eNTU(3000,0.99,:counterflow) ≈ 1.0 atol = 1e-5
end