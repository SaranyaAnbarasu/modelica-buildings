within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model BoilerWithEconomizerDynamic
  extends Modelica.Icons.Example;

  //data file
  parameter String data = ("modelica://Buildings/Resources/Data/GEDCalibration/BoilerWithEconCase23328000.mos");

  // Medium declarations
  package MediumWat =
      Buildings.Media.Specialized.Water.TemperatureDependentDensity (
      p_default=101325,
      T_default=273.15+100)
        "Water medium - port_a (inlet)";
  package MediumSte = Buildings.Media.Steam (
    p_default=300000,
    T_default=273.15+200,
    h_default=2700000)
     "Steam medium - port_b (oulet)";
  package MediumAir = Buildings.Media.CombustionAir
     "Combustion air medium";

  // Boiler nominal conditions
  parameter Modelica.Units.SI.AbsolutePressure p_nominal=917003
    "Nominal pressure";
  parameter Modelica.Units.SI.Temperature T_nominal=
      MediumSte.saturationTemperature(p_nominal)
    "Nominal saturation temperature";
  parameter Modelica.Units.SI.Power Q_flow_nominal=17496340 "Nominal power";
  parameter Modelica.Units.SI.SpecificEnthalpy dh_nominal=
      MediumSte.specificEnthalpy(MediumSte.setState_pTX(
      p=p_nominal,
      T=T_nominal,
      X=MediumSte.X_default)) "Nominal change in enthalpy";
  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=Q_flow_nominal/
      dh_nominal/2 "Nominal mass flow rate";
  parameter Modelica.Units.SI.PressureDifference dp_nominal=800000
    "Pressure drop at m_flow_nominal";

  parameter Modelica.Units.SI.Temperature T_exh_nominal(displayUnit="K") = 423
    "Exhaust temperature used to compute nominal efficiency";
  parameter Real FA_ratio=1.15
    "Fuel air ratio, alpha (20% excess air, FA_ratio = 1.20)";
  //parameter Fluid.Data.Fuels.Generic fue=fue(
   //   h=46402971,
    //  d=800,
    //  mCO2=2.2) "Fuel type";
  parameter Modelica.Units.SI.Volume V=50 "Total internal volume of boiler";
  parameter Modelica.Units.SI.Volume V_com=20
    "Total internal volume of combustion";
  parameter Modelica.Units.SI.Mass mDry=1.5E-3*Q_flow_nominal
    "Mass of boiler that will be lumped to water heat capacity";
  parameter Modelica.Media.Interfaces.Types.AbsolutePressure p_start=1000000
    "Start value of pressure";

  //economizer
  parameter Real r_nominal=2/3
    "Ratio between air-side and water-side convective heat transfer (hA-value) at nominal condition";
    parameter Modelica.Units.SI.Power Q_eco_flow_nominal=1324000 "Nominal power of economizer";

  parameter Modelica.Units.SI.Temperature T_a2_nominal(displayUnit="K") = 377
    "Nominal temperature at port a2";
  parameter Modelica.Units.SI.Temperature T_a1_nominal(displayUnit="K") = 549
    "Nominal temperature at port a1";

  parameter Modelica.Units.SI.MassFlowRate m1_flow_nominal=m_flow_nominal
    "Nominal mass flow rate between port a1 and b1";

  parameter Modelica.Units.SI.MassFlowRate m2_flow_nominal=6.36
    "Nominal mass flow rate between port a2 and b2";

  parameter Modelica.Units.SI.PressureDifference dp1_nominal=20000
    "Pressure drop at m_flow_nominal between port a1 and b1";

  parameter Modelica.Units.SI.PressureDifference dp2_nominal=39989.6
    "Pressure drop at m_flow_nominal between port a2 and b2";

  parameter Modelica.Units.SI.PressureDifference dpValve_nominal=50000
    "Nominal pressure drop of fully open valve, used if CvData=Buildings.Fluid.Types.CvTypes.OpPoint";

// Dynamics
  parameter Modelica.Fluid.Types.Dynamics energyDynamics=
    Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
    "Type of energy balance: dynamic (3 initialization options) or steady state"
    annotation (Evaluate=true,Dialog(tab="Dynamics",group="Equations"));
  parameter Modelica.Fluid.Types.Dynamics massDynamics=energyDynamics
    "Type of mass balance: dynamic (3 initialization options) or steady state"
    annotation (Evaluate=true,Dialog(tab="Dynamics",group="Equations"));

  Fluid.Movers.FlowControlled_m_flow fwPum(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    addPowerToMedium=false,
    nominalValuesDefineDefaultPressureCurve=true,
    dp_nominal=dp_nominal) "Feed water pump"
    annotation (Placement(transformation(extent={{-70,10},{-50,30}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    plr(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={4},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-140,-80},{-120,-60}})));
  Modelica.Blocks.Sources.Constant Tair_in(k=294.15)
    annotation (Placement(transformation(extent={{-392,4},{-372,24}})));
  Fluid.Sources.Boundary_pT pro(
    redeclare package Medium = MediumAir,
    use_T_in=true,
    T=573.15,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-100,-110},{-80,-90}})));
  HeatTransfer.Sources.FixedTemperature TAmb(T(displayUnit="K") = 294)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{120,-80},{100,-60}})));
  Fluid.Sources.Boundary_pT exh(
    redeclare package Medium = MediumAir,
    p=100000,
    T(displayUnit="K") = 560,
    nPorts=1) "Source" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={66,70})));

  Fluid.Sources.Boundary_pT sin(
    redeclare package Medium = MediumSte,
    p(displayUnit="bar") = 900000,
    T=453.15,
    nPorts=1) "Sink"
    annotation (Placement(transformation(extent={{154,-118},{134,-98}})));
  Modelica.Blocks.Sources.Constant ecoValSig(k=1)
    annotation (Placement(transformation(extent={{-140,70},{-120,90}})));
  Fluid.Sensors.TemperatureTwoPort senTem(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={66,38})));
  Fluid.Sensors.TemperatureTwoPort senTem1(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={66,-52})));
  Fluid.Sensors.TemperatureTwoPort senTem2(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={72,-108})));
  Fluid.Sources.Boundary_pT sou(
    redeclare package Medium = MediumWat,
    p=100000,
    use_T_in=true,
    T=303.15,
    nPorts=1) "Source" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-90,20})));
  Modelica.Blocks.Sources.Constant Tfw_in(k=377) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-344,-110})));
  Components.BoilerPolynomialFurnaceHeatBalance boi(
    m1_flow_nominal=m1_flow_nominal,
    m2_flow_nominal=m2_flow_nominal,
    show_T=true,
    energyDynamics=energyDynamics,
    massDynamics=massDynamics,
    p_start=p_start,
    redeclare package MediumAir = MediumAir,
    redeclare package MediumWat = MediumWat,
    redeclare package MediumSte = MediumSte,
    mDry=mDry,
    m_flow_nominal=10,
    dp_nominal=400000,
    Q_flow_nominal=Q_flow_nominal,
    fue=Buildings.Fluid.Data.Fuels.NaturalGasLowerHeatingValue(h=46402971),
    UA=0.01*Q_flow_nominal/100,
    V=V,
    V_com=V_com,
    FA_ratio=FA_ratio,
    T_exh_nominal(displayUnit="K") = T_exh_nominal)
    annotation (Placement(transformation(extent={{12,-114},{32,-94}})));
  Fluid.HeatExchangers.DryCoilEffectivenessNTU eco(
    redeclare package Medium1 = MediumAir,
    redeclare package Medium2 = MediumWat,
    allowFlowReversal1=false,
    allowFlowReversal2=false,
    m1_flow_nominal=m1_flow_nominal,
    m2_flow_nominal=m2_flow_nominal,
    show_T=true,
    dp1_nominal=dp1_nominal,
    dp2_nominal=dp2_nominal,
    configuration=Buildings.Fluid.Types.HeatExchangerConfiguration.CounterFlow,
    use_Q_flow_nominal=false,
    T_a1_nominal(displayUnit="K"),
    T_a2_nominal(displayUnit="K"),
    eps_nominal=0.85,
    r_nominal=0.5) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={60,-10})));

  Fluid.Actuators.Valves.TwoWayLinear val(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{6,10},{26,30}})));
  Fluid.Actuators.Valves.TwoWayLinear val1(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{26,-50},{6,-30}})));
  Fluid.Actuators.Valves.TwoWayLinear val2(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
                           annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-26,0})));
  Modelica.Blocks.Sources.Constant valSig(k=0)
    annotation (Placement(transformation(extent={{-140,-30},{-120,-10}})));
  Modelica.Blocks.Math.Gain m_fw_nom(k=0.969)
    annotation (Placement(transformation(extent={{-100,42},{-80,62}})));

  Modelica.Blocks.Sources.Sine inp(
    amplitude=-0.5,
    f=1/86400,
    phase=3.1415926535898,
    offset=0.5)
    "Input signal"
    annotation (Placement(transformation(extent={{-396,62},{-376,82}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    T_flu(
    tableOnFile=true,
    table=fill(
        0.0,
        0,
        2),
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={5},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    "Combustion air and flue gas mix temperature"
    annotation (Placement(transformation(extent={{-140,-106},{-120,-86}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    mFloFW1(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={2},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-140,42},{-120,62}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    outputs(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns=2:12,
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-354,-82},{-334,-62}})));
  Modelica.Blocks.Sources.CombiTimeTable T_FW(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={5},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-140,6},{-120,26}})));
  Fluid.FixedResistances.PressureDrop res1(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m1_flow_nominal,
    dp_nominal=10000)
    annotation (Placement(transformation(extent={{98,-118},{118,-98}})));
equation

  connect(eco.port_a2, val.port_b) annotation (Line(
      points={{54,0},{54,20},{26,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_a, eco.port_b2) annotation (Line(
      points={{26,-40},{54,-40},{54,-20}},
      color={0,127,255},
      thickness=0.5));
  connect(val2.port_a, fwPum.port_b) annotation (Line(
      points={{-26,10},{-26,20},{-50,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val2.port_b, boi.port_a2) annotation (Line(points={{-26,-10},{-26,-108},
          {12,-108}}, color={0,127,255},
      thickness=0.5));
  connect(val.port_a, fwPum.port_b) annotation (Line(
      points={{6,20},{-50,20}},
      color={0,127,255},
      thickness=0.5));
  connect(sou.ports[1], fwPum.port_a) annotation (Line(
      points={{-80,20},{-70,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_b, boi.port_a2) annotation (Line(
      points={{6,-40},{-26,-40},{-26,-108},{12,-108}},
      color={0,127,255},
      thickness=0.5));
  connect(TAmb.port, boi.heatPort)
    annotation (Line(points={{100,-70},{22,-70},{22,-94}},color={191,0,0}));
  connect(boi.port_b1, senTem1.port_a) annotation (Line(
      points={{32,-100},{66,-100},{66,-62}},
      color={0,0,0},
      thickness=0.5));
  connect(senTem1.port_b, eco.port_a1) annotation (Line(
      points={{66,-42},{66,-20}},
      color={0,0,0},
      thickness=0.5));
  connect(pro.ports[1], boi.port_a1)
    annotation (Line(points={{-80,-100},{12,-100}},
                                                  color={148,145,145},
      thickness=0.5));
  connect(boi.port_b2, senTem2.port_a) annotation (Line(
      points={{32,-108},{62,-108}},
      color={238,46,47},
      thickness=0.5));
  connect(exh.ports[1], senTem.port_b) annotation (Line(
      points={{66,60},{66,48}},
      color={0,0,0},
      thickness=0.5));
  connect(senTem.port_a, eco.port_b1) annotation (Line(
      points={{66,28},{66,0}},
      color={0,0,0},
      thickness=0.5));
  connect(ecoValSig.y, val.y) annotation (Line(
      points={{-119,80},{16,80},{16,32}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig.y, val1.y) annotation (Line(
      points={{-119,80},{0,80},{0,4},{16,4},{16,-28}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(valSig.y, val2.y) annotation (Line(
      points={{-119,-20},{-40,-20},{-40,6.66134e-16},{-38,6.66134e-16}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(m_fw_nom.y, fwPum.m_flow_in)
    annotation (Line(points={{-79,52},{-60,52},{-60,32}},  color={0,0,127},
      pattern=LinePattern.Dash));
  connect(T_flu.y[1], pro.T_in)
    annotation (Line(points={{-119,-96},{-102,-96}},color={0,0,127},
      pattern=LinePattern.Dash));
  connect(plr.y[1], boi.y) annotation (Line(points={{-119,-70},{0,-70},{0,-96},{
          10,-96}},  color={0,0,127},
      pattern=LinePattern.Dash));
  connect(mFloFW1.y[1], m_fw_nom.u)
    annotation (Line(points={{-119,52},{-102,52}}, color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senTem2.port_b, res1.port_a) annotation (Line(
      points={{82,-108},{98,-108}},
      color={238,46,47},
      thickness=0.5));
  connect(res1.port_b, sin.ports[1]) annotation (Line(
      points={{118,-108},{134,-108}},
      color={238,46,47},
      thickness=0.5));
  connect(T_FW.y[1], sou.T_in) annotation (Line(
      points={{-119,16},{-102,16}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  annotation (
    __Dymola_Commands(file="modelica://Buildings/Resources/Data/GEDCalibration/Scripts/BoilerWithEconomizerDynamic.mos"
        "Simulate and plot"),
    experiment(
      StartTime=23328000,
      StopTime=25056000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,
            100}})));
end BoilerWithEconomizerDynamic;
