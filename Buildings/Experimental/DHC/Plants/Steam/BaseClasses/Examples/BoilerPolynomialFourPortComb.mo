within Buildings.Experimental.DHC.Plants.Steam.BaseClasses.Examples;
model BoilerPolynomialFourPortComb
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
  parameter Modelica.Units.SI.AbsolutePressure p_nominal = 917003
    "Nominal pressure";
  parameter Modelica.Units.SI.Temperature T_nominal=
    MediumSte.saturationTemperature(p_nominal)
    "Nominal saturation temperature";
  parameter Modelica.Units.SI.Power Q_flow_nominal = 29306000 "Nominal power";
  parameter Modelica.Units.SI.SpecificEnthalpy dh_nominal=
    MediumSte.specificEnthalpy(
      MediumSte.setState_pTX(p=p_nominal, T=T_nominal, X=MediumSte.X_default))
    "Nominal change in enthalpy";
  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=
    Q_flow_nominal/dh_nominal/2
    "Nominal mass flow rate";
  parameter Modelica.Units.SI.PressureDifference dp_nominal = 800000
    "Pressure drop at m_flow_nominal";

  Buildings.Fluid.Sources.Boundary_pT sin(
    redeclare package Medium = MediumSte,
    p(displayUnit="bar") = 900000,
    T=453.15,
    nPorts=1)
    "Sink"
    annotation (Placement(transformation(extent={{80,10},{60,30}})));
  Buildings.Fluid.Sources.Boundary_pT sou(
    redeclare package Medium = MediumWat,
    T=303.15,
    nPorts=1)
    "Source"
    annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
  Buildings.HeatTransfer.Sources.FixedTemperature TAmb(T=303.15)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{-30,102},{-10,122}})));
  BoilerPolynomialExhaustCombustiongas
                          boiDyn(
    m1_flow_nominal=1,
    m2_flow_nominal=m_flow_nominal,
    redeclare package MediumAir = MediumAir,
    redeclare package MediumSte = MediumSte,
    redeclare package MediumWat = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    Q_flow_nominal=Q_flow_nominal,
    a={0.8},
    fue=Buildings.Fluid.Data.Fuels.NaturalGasHigherHeatingValue(),
    dp_nominal=dp_nominal,
    UA=0.05*Q_flow_nominal/100,
    V=25,
    V_com=21.94)           "Steam boiler with dynamic balance"
    annotation (Placement(transformation(extent={{-10,42},{10,62}})));
  Fluid.Sources.Boundary_pT           sou1(
    redeclare package Medium = MediumAir,
    use_T_in=true,
    T=573.15,
    nPorts=1)
    "Source"
    annotation (Placement(transformation(extent={{-80,70},{-60,90}})));
  Fluid.Sources.Boundary_pT           sou2(redeclare package Medium = MediumAir,
      nPorts=1)
    "Source"
    annotation (Placement(transformation(extent={{80,70},{60,90}})));
  Modelica.Blocks.Sources.Constant
                               const(k=273.15 + 48)
    annotation (Placement(transformation(extent={{-120,74},{-100,94}})));
  Fluid.Movers.FlowControlled_m_flow fwPum(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m_flow_nominal,
    addPowerToMedium=false,
    nominalValuesDefineDefaultPressureCurve=true,
    dp_nominal=dp_nominal) "Feed water pump"
    annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
  Modelica.Blocks.Sources.Constant
                               const1(k=4.8)
    annotation (Placement(transformation(extent={{-120,30},{-100,50}})));
  Modelica.Blocks.Sources.Pulse pulse(
    amplitude=1,
    width=30,
    period=500)
    annotation (Placement(transformation(extent={{-120,114},{-100,134}})));
  Modelica.Blocks.Sources.RealExpression QFue(y=boiDyn.QFue_flow)
    "Fuel heat flow rate"
    annotation (Placement(transformation(extent={{-170,-66},{-150,-46}})));
  Modelica.Blocks.Sources.RealExpression QWat(y=boiDyn.QWat_flow)
    "Water heat flow rate"
    annotation (Placement(transformation(extent={{-170,-116},{-150,-96}})));
  Modelica.Blocks.Sources.RealExpression QLoss(y=boiDyn.heatPort.Q_flow)
    "Loss from boiler casing"
    annotation (Placement(transformation(extent={{-62,-118},{-42,-98}})));
  Modelica.Blocks.Sources.RealExpression QFlue(y=boiDyn.Q_flow_exh.y)
    "Heat losses into flue gas"
    annotation (Placement(transformation(extent={{-62,-68},{-42,-48}})));
  Modelica.Blocks.Sources.RealExpression Qcomb(y=((boiDyn.port_b1.h_outflow -
        boiDyn.port_a1.h_outflow)*boiDyn.port_a1.m_flow))
    "Q calcuated based on the mflow and enthalpy of the combustion side ports"
    annotation (Placement(transformation(extent={{66,-68},{86,-48}})));
  Modelica.Blocks.Sources.RealExpression Qevap(y=((boiDyn.port_b2.h_outflow -
        boiDyn.port_a2.h_outflow)*boiDyn.port_b2.m_flow))
    "Q calcuated based on the mflow and enthalpy of the combustion side ports"
    annotation (Placement(transformation(extent={{66,-118},{86,-98}})));
  Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
    annotation (Placement(transformation(extent={{-132,-66},{-112,-46}})));
  Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
    annotation (Placement(transformation(extent={{-132,-116},{-112,-96}})));
  Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
    annotation (Placement(transformation(extent={{-20,-68},{0,-48}})));
  Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
    annotation (Placement(transformation(extent={{-20,-118},{0,-98}})));
  Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
    annotation (Placement(transformation(extent={{114,-68},{134,-48}})));
  Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
    annotation (Placement(transformation(extent={{114,-118},{134,-98}})));
  Modelica.Blocks.Math.Division perQWat
    annotation (Placement(transformation(extent={{-96,-90},{-76,-110}})));
  Modelica.Blocks.Math.Division perQLoss
    annotation (Placement(transformation(extent={{14,-92},{34,-112}})));
  Modelica.Blocks.Math.Division perQFlue
    annotation (Placement(transformation(extent={{16,-74},{36,-54}})));
  Modelica.Blocks.Math.Division perQcomb
    annotation (Placement(transformation(extent={{150,-74},{170,-54}})));
  Modelica.Blocks.Math.Division perQeva
    annotation (Placement(transformation(extent={{150,-92},{170,-112}})));
  Modelica.Blocks.Math.Add diffQevap(k1=-1)
    annotation (Placement(transformation(extent={{32,-146},{12,-166}})));
  Modelica.Blocks.Math.Add diffQcomb(k1=-1)
    annotation (Placement(transformation(extent={{90,-146},{70,-166}})));
equation
  connect(TAmb.port, boiDyn.heatPort)
    annotation (Line(points={{-10,112},{0,112},{0,62}},   color={191,0,0}));
  connect(boiDyn.port_b2, sin.ports[1]) annotation (Line(points={{10,48},{20,48},
          {20,20},{60,20}},   color={0,127,255}));
  connect(sou1.ports[1], boiDyn.port_a1) annotation (Line(points={{-60,80},{-42,
          80},{-42,56},{-10,56}},
                              color={0,127,255}));
  connect(boiDyn.port_b1, sou2.ports[1]) annotation (Line(points={{10,56},{20,
          56},{20,80},{60,80}},
                            color={0,127,255}));
  connect(const.y, sou1.T_in)
    annotation (Line(points={{-99,84},{-82,84}}, color={0,0,127}));
  connect(sou.ports[1], fwPum.port_a)
    annotation (Line(points={{-60,20},{-40,20}},   color={0,127,255}));
  connect(fwPum.port_b, boiDyn.port_a2) annotation (Line(points={{-20,20},{-16,
          20},{-16,48},{-10,48}},  color={0,127,255}));
  connect(const1.y, fwPum.m_flow_in)
    annotation (Line(points={{-99,40},{-30,40},{-30,32}},    color={0,0,127}));
  connect(pulse.y, boiDyn.y) annotation (Line(points={{-99,124},{-40,124},{-40,
          60},{-12,60}}, color={0,0,127}));
  connect(QFue.y, iQfue.u)
    annotation (Line(points={{-149,-56},{-134,-56}}, color={0,0,127}));
  connect(QWat.y, iQWat.u)
    annotation (Line(points={{-149,-106},{-134,-106}}, color={0,0,127}));
  connect(QFlue.y, iQFlue.u)
    annotation (Line(points={{-41,-58},{-22,-58}}, color={0,0,127}));
  connect(QLoss.y, iQLoss.u)
    annotation (Line(points={{-41,-108},{-22,-108}}, color={0,0,127}));
  connect(iQComb.u, Qcomb.y)
    annotation (Line(points={{112,-58},{87,-58}}, color={0,0,127}));
  connect(iQevap.u, Qevap.y)
    annotation (Line(points={{112,-108},{87,-108}}, color={0,0,127}));
  connect(iQWat.y, perQWat.u1)
    annotation (Line(points={{-111,-106},{-98,-106}}, color={0,0,127}));
  connect(iQfue.y, perQWat.u2) annotation (Line(
      points={{-111,-56},{-104,-56},{-104,-94},{-98,-94}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(iQLoss.y, perQLoss.u1)
    annotation (Line(points={{1,-108},{12,-108}}, color={0,0,127}));
  connect(iQfue.y, perQLoss.u2) annotation (Line(
      points={{-111,-56},{-104,-56},{-104,-82},{6,-82},{6,-96},{12,-96}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(iQFlue.y, perQFlue.u1)
    annotation (Line(points={{1,-58},{14,-58}}, color={0,0,127}));
  connect(iQfue.y, perQFlue.u2) annotation (Line(
      points={{-111,-56},{-104,-56},{-104,-82},{6,-82},{6,-70},{14,-70}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(iQComb.y, perQcomb.u1)
    annotation (Line(points={{135,-58},{148,-58}}, color={0,0,127}));
  connect(iQfue.y, perQcomb.u2) annotation (Line(
      points={{-111,-56},{-104,-56},{-104,-82},{140,-82},{140,-70},{148,-70}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(iQevap.y, perQeva.u1)
    annotation (Line(points={{135,-108},{148,-108}}, color={0,0,127}));
  connect(iQfue.y, perQeva.u2) annotation (Line(
      points={{-111,-56},{-104,-56},{-104,-82},{140,-82},{140,-96},{148,-96}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(perQcomb.y, diffQcomb.u2) annotation (Line(points={{171,-64},{178,-64},
          {178,-150},{92,-150}}, color={0,0,127}));
  connect(perQeva.y, diffQevap.u2) annotation (Line(points={{171,-102},{170,-102},
          {170,-142},{34,-142},{34,-150}}, color={0,0,127}));
  connect(perQWat.y, diffQevap.u1) annotation (Line(points={{-75,-100},{-68,-100},
          {-68,-176},{46,-176},{46,-162},{34,-162}}, color={244,125,35}));
  connect(perQFlue.y, diffQcomb.u1) annotation (Line(points={{37,-64},{60,-64},{
          60,-176},{102,-176},{102,-162},{92,-162}}, color={244,125,35}));
  annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Plants/Steam/BaseClasses/Examples/BoilerPolynomialFourPort.mos"
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
</html>"),
    Diagram(coordinateSystem(extent={{-180,-140},{180,140}}), graphics={
          Rectangle(
          extent={{-180,-40},{180,-180}},
          lineColor={28,108,200},
          fillColor={174,179,179},
          fillPattern=FillPattern.Solid), Text(
          extent={{-176,-26},{-70,-46}},
          textColor={28,108,200},
          textString="Steady state verification (Heat balance)")}),
    Icon(coordinateSystem(extent={{-180,-140},{180,140}})));
end BoilerPolynomialFourPortComb;
