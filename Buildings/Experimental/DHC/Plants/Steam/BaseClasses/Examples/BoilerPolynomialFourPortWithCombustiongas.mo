within Buildings.Experimental.DHC.Plants.Steam.BaseClasses.Examples;
model BoilerPolynomialFourPortWithCombustiongas
  "Example model for the steam boiler with a polynomial efficiency curve"
  extends Modelica.Icons.Example;

  // Medium declarations
  package MediumWat =
      Buildings.Media.Specialized.Water.TemperatureDependentDensity
    "Water medium - port_a (inlet)";
  package MediumSte = Buildings.Media.Steam
     "Steam medium - port_b (oulet)";
  package MediumAir = Buildings.Media.CombustionAir
     "Combustion air medium";

  // Nominal conditions
  parameter Modelica.Units.SI.AbsolutePressure p_nominal = 300000
    "Nominal pressure";
  parameter Modelica.Units.SI.Temperature T_nominal=
    MediumSte.saturationTemperature(p_nominal)
    "Nominal saturation temperature";
  parameter Modelica.Units.SI.Power Q_flow_nominal = 50000 "Nominal power";
  parameter Modelica.Units.SI.SpecificEnthalpy dh_nominal=
    MediumSte.specificEnthalpy(
      MediumSte.setState_pTX(p=p_nominal, T=T_nominal, X=MediumSte.X_default))
    "Nominal change in enthalpy";
  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=
    Q_flow_nominal/dh_nominal/2
    "Nominal mass flow rate";
  parameter Modelica.Units.SI.PressureDifference dp_nominal = 3000
    "Pressure drop at m_flow_nominal";

  Modelica.Blocks.Sources.TimeTable y(table=[0,0; 1200,1; 1200,0; 2000,0; 2000,
        1; 3600,0])
    "Load ratio"
    annotation (Placement(transformation(extent={{-70,48},{-50,68}})));
  Buildings.Fluid.Sources.Boundary_pT sin(
    redeclare package Medium = MediumSte,
    p(displayUnit="bar") = 300000,
    T=423.15,
    nPorts=1)
    "Sink"
    annotation (Placement(transformation(extent={{80,-40},{60,-20}})));
  Buildings.Fluid.Sources.Boundary_pT sou(
    redeclare package Medium = MediumWat,
    p=300000 + dp_nominal,
    T=303.15,
    nPorts=1)
    "Source"
    annotation (Placement(transformation(extent={{-80,-42},{-60,-22}})));
  Buildings.HeatTransfer.Sources.FixedTemperature TAmb(T=288.15)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{-30,48},{-10,68}})));
  BoilerPolynomialExhaustCombustiongas
                          boiDyn(
    m1_flow_nominal=0.2,
    m2_flow_nominal=m_flow_nominal,
    redeclare package MediumAir = MediumAir,
    redeclare package MediumSte = MediumSte,
    redeclare package MediumWat = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    Q_flow_nominal=Q_flow_nominal,
    fue=Buildings.Fluid.Data.Fuels.NaturalGasLowerHeatingValue(),
    dp_nominal=dp_nominal) "Steam boiler with dynamic balance"
    annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
  Fluid.Sources.Boundary_pT           sou1(
    redeclare package Medium = MediumAir,
    T=333.15,
    nPorts=1)
    "Source"
    annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
  Fluid.Sources.Boundary_pT           sou2(redeclare package Medium = MediumAir,
      nPorts=1)
    "Source"
    annotation (Placement(transformation(extent={{80,10},{60,30}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=273.15 + 20,
    duration=1000,
    startTime=500)
    annotation (Placement(transformation(extent={{-120,14},{-100,34}})));
equation
  connect(TAmb.port, boiDyn.heatPort)
    annotation (Line(points={{-10,58},{0,58},{0,8}},      color={191,0,0}));
  connect(y.y, boiDyn.y)
    annotation (Line(points={{-49,58},{-40,58},{-40,6},{-12,6}},
                   color={0,0,127}));
  connect(boiDyn.port_b2, sin.ports[1]) annotation (Line(points={{10,-6},{54,-6},
          {54,-30},{60,-30}}, color={0,127,255}));
  connect(sou1.ports[1], boiDyn.port_a1) annotation (Line(points={{-60,20},{-42,
          20},{-42,2},{-10,2}},
                              color={0,127,255}));
  connect(boiDyn.port_b1, sou2.ports[1]) annotation (Line(points={{10,2},{54,2},
          {54,20},{60,20}}, color={0,127,255}));
  connect(sou.ports[1], boiDyn.port_a2) annotation (Line(points={{-60,-32},{-16,
          -32},{-16,-6},{-10,-6}}, color={0,127,255}));
  annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Plants/Steam/BaseClasses/Examples/BoilerPolynomial.mos"
        "Simulate and plot"),
    experiment(Tolerance=1e-6, StopTime=3600),
    Documentation(info="<html>
<p>
This example demonstrates the open loop response of the 
steam boiler model. The dynamic boiler includes a control 
signal that is first a ramp from <i>0</i> to <i>1</i>, 
followed by a step that switches the boiler off and then 
on again. The steady boiler is only dependent on the fluid
flow.
</p>
</html>", revisions="<html>
<ul>
<li>
July 23, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"));
end BoilerPolynomialFourPortWithCombustiongas;
