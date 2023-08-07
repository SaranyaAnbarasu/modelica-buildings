within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model BoilerWithEconomizerSteadyState
  extends Modelica.Icons.Example;

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
  parameter Modelica.Units.SI.Volume V_com=10
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
    annotation (Placement(transformation(extent={{-86,8},{-66,28}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    plr(table=[0.0,0; 900,0.25; 1800,0.5; 2700,0.75;
        3600,1])
    annotation (Placement(transformation(extent={{-166,-64},{-146,-44}})));
  Modelica.Blocks.Sources.Constant Tair_in(k=294.15)
    annotation (Placement(transformation(extent={{-392,4},{-372,24}})));
  Fluid.Sources.Boundary_pT pro(
    redeclare package Medium = MediumAir,
    use_T_in=true,
    T=573.15,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-86,-104},{-66,-84}})));
  HeatTransfer.Sources.FixedTemperature TAmb(T(displayUnit="K") = 294)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{54,-72},{34,-52}})));
  Fluid.Sources.Boundary_pT exh(
    redeclare package Medium = MediumAir,
    p=100000,
    T(displayUnit="K") = 560,
    nPorts=1) "Source" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={66,68})));

  Fluid.Sources.Boundary_pT sin(
    redeclare package Medium = MediumSte,
    p(displayUnit="bar") = 900000,
    T=453.15,
    nPorts=1) "Sink"
    annotation (Placement(transformation(extent={{154,-112},{134,-92}})));
  Modelica.Blocks.Sources.Constant ecoValSig(k=1)
    annotation (Placement(transformation(extent={{-166,68},{-146,88}})));
  Fluid.Sensors.TemperatureTwoPort senTem(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={66,36})));
  Fluid.Sensors.TemperatureTwoPort senTem1(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={66,-32})));
  Fluid.Sensors.TemperatureTwoPort senTem2(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={104,-102})));
  Fluid.Sources.Boundary_pT sou(
    redeclare package Medium = MediumWat,
    p=100000,
    use_T_in=true,
    T=303.15,
    nPorts=1) "Source" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-116,18})));
  Modelica.Blocks.Sources.Constant Tfw_in(k=377) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-156,14})));
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
    fue(
      h=42858187,
      d=800000,
      mCO2=2.2),
    UA=0.01*Q_flow_nominal/100,
    V=V,
    V_com=V_com,
    FA_ratio=FA_ratio,
    T_exh_nominal(displayUnit="K") = T_exh_nominal)
    annotation (Placement(transformation(extent={{12,-108},{32,-88}})));
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
    r_nominal=0.8) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={60,2})));

  Fluid.Actuators.Valves.TwoWayLinear val(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{6,8},{26,28}})));
  Fluid.Actuators.Valves.TwoWayLinear val1(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{26,-28},{6,-8}})));
  Fluid.Actuators.Valves.TwoWayLinear val2(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
                           annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-26,0})));
  Modelica.Blocks.Sources.Constant valSig(k=0)
    annotation (Placement(transformation(extent={{-166,-32},{-146,-12}})));
  Modelica.Blocks.Math.Gain m_fw_nom(k=6.30)
    annotation (Placement(transformation(extent={{-120,40},{-100,60}})));

  Modelica.Blocks.Sources.Sine inp(
    amplitude=-0.5,
    f=1/86400,
    phase=3.1415926535898,
    offset=0.5)
    "Input signal"
    annotation (Placement(transformation(extent={{-396,62},{-376,82}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    T_flu(table=[0.0,300; 900,302; 1800,
        304; 2700,306; 3600,308])
                          "Combustion air and flue gas mix temperature"
    annotation (Placement(transformation(extent={{-166,-100},{-146,-80}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    mFloFW1(table=[0.0,0; 900,0.25; 1800,0.5; 2700,
        0.75; 3600,1])
    annotation (Placement(transformation(extent={{-166,40},{-146,60}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    outputs(table=fill(
              0.0,
              0,
              2))
    annotation (Placement(transformation(extent={{-260,-210},{-240,-190}})));
equation

  connect(eco.port_a2, val.port_b) annotation (Line(
      points={{54,12},{54,18},{26,18}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_a, eco.port_b2) annotation (Line(
      points={{26,-18},{54,-18},{54,-8}},
      color={0,127,255},
      thickness=0.5));
  connect(val2.port_a, fwPum.port_b) annotation (Line(
      points={{-26,10},{-26,18},{-66,18}},
      color={0,127,255},
      thickness=0.5));
  connect(val2.port_b, boi.port_a2) annotation (Line(points={{-26,-10},{-26,-102},
          {12,-102}}, color={0,127,255}));
  connect(val.port_a, fwPum.port_b) annotation (Line(
      points={{6,18},{-66,18}},
      color={0,127,255},
      thickness=0.5));
  connect(sou.ports[1], fwPum.port_a) annotation (Line(
      points={{-106,18},{-86,18}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_b, boi.port_a2) annotation (Line(
      points={{6,-18},{-26,-18},{-26,-102},{12,-102}},
      color={0,127,255},
      thickness=0.5));
  connect(TAmb.port, boi.heatPort)
    annotation (Line(points={{34,-62},{22,-62},{22,-88}}, color={191,0,0}));
  connect(boi.port_b1, senTem1.port_a) annotation (Line(
      points={{32,-94},{66,-94},{66,-42}},
      color={0,0,0},
      thickness=0.5));
  connect(senTem1.port_b, eco.port_a1) annotation (Line(
      points={{66,-22},{66,-8}},
      color={0,0,0},
      thickness=0.5));
  connect(sou.T_in, Tfw_in.y)
    annotation (Line(points={{-128,14},{-145,14}}, color={0,0,127}));
  connect(pro.ports[1], boi.port_a1)
    annotation (Line(points={{-66,-94},{12,-94}}, color={0,127,255}));
  connect(boi.port_b2, senTem2.port_a) annotation (Line(
      points={{32,-102},{94,-102}},
      color={238,46,47},
      thickness=0.5));
  connect(senTem2.port_b, sin.ports[1]) annotation (Line(
      points={{114,-102},{134,-102}},
      color={238,46,47},
      thickness=0.5));
  connect(exh.ports[1], senTem.port_b) annotation (Line(
      points={{66,58},{66,46}},
      color={0,0,0},
      thickness=0.5));
  connect(senTem.port_a, eco.port_b1) annotation (Line(
      points={{66,26},{66,12}},
      color={0,0,0},
      thickness=0.5));
  connect(ecoValSig.y, val.y) annotation (Line(
      points={{-145,78},{16,78},{16,30}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig.y, val1.y) annotation (Line(
      points={{-145,78},{0,78},{0,2},{16,2},{16,-6}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(valSig.y, val2.y) annotation (Line(
      points={{-145,-22},{-44,-22},{-44,0},{-38,0}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(plr.y, boi.y) annotation (Line(
      points={{-145,-54},{-10,-54},{-10,-90},{10,-90}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(T_flu.y, pro.T_in)
    annotation (Line(points={{-145,-90},{-88,-90}}, color={0,0,127}));
  connect(mFloFW1.y, m_fw_nom.u)
    annotation (Line(points={{-145,50},{-122,50}}, color={0,0,127}));
  connect(m_fw_nom.y, fwPum.m_flow_in)
    annotation (Line(points={{-99,50},{-76,50},{-76,30}},  color={0,0,127}));
  annotation (
    __Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Components/HeatBalanceBoiler.mos"
        "Simulate and plot"),
    experiment(Tolerance=1e-6, StopTime=3600),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,
            100}})));
end BoilerWithEconomizerSteadyState;
