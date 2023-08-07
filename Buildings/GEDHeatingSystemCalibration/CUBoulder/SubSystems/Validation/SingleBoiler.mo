within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model SingleBoiler "Example model to demonstrate the single-boiler steam plant 
  in a single closed loop"
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

  // Nominal conditions
  parameter Modelica.Units.SI.AbsolutePressure p_nominal = 917003
    "Nominal pressure";
  parameter Modelica.Units.SI.Temperature T_nominal=
    MediumSte.saturationTemperature(p_nominal)
    "Nominal saturation temperature";
  parameter Modelica.Units.SI.Power Q_flow_nominal = 17496340 "Nominal power";
  parameter Modelica.Units.SI.SpecificEnthalpy dh_nominal=
    MediumSte.specificEnthalpy(
      MediumSte.setState_pTX(p=p_nominal, T=T_nominal, X=MediumSte.X_default))
    "Nominal change in enthalpy";
  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=
    Q_flow_nominal/dh_nominal/2
    "Nominal mass flow rate";
  parameter Modelica.Units.SI.PressureDifference dp_nominal = 800000
    "Pressure drop at m_flow_nominal";

  parameter Modelica.Units.SI.AbsolutePressure pSat=800000
    "Saturation pressure";
  parameter Modelica.Units.SI.Temperature TSat=
     MediumSte.saturationTemperature(pSat)
     "Saturation temperature";

  parameter Modelica.Units.SI.PressureDifference dpPip=6000
    "Pressure drop in the condensate return pipe";
  // pumps
  parameter Buildings.Fluid.Movers.Data.Generic perPumFW(
    pressure(
      V_flow=m_flow_nominal*1000*{0.4,0.6,0.8,1.0},
      dp=(pSat-101325)*{1.34,1.27,1.17,1.0}))
    "Performance data for feedwater pump";
    parameter Buildings.Fluid.Movers.Data.Generic perPumCNR(
   pressure(
     V_flow=m_flow_nominal*1000*{0,1,2},
     dp=dpPip*{2,1,0}))
    "Performance data for condensate return pumps";

      // Boiler controller
  parameter Modelica.Blocks.Types.SimpleController controllerTypeBoi=
    Modelica.Blocks.Types.SimpleController.PI "Type of controller"
    annotation (Dialog(tab="Control", group="Boiler"));
  parameter Real kBoi(min=0) = 10 "Gain of controller"
    annotation (Dialog(tab="Control", group="Boiler"));
  parameter Modelica.Units.SI.Time TiBoi(min=Modelica.Constants.small)=120
    "Time constant of Integrator block"
     annotation (Dialog(enable=
          controllerTypeBoi == Modelica.Blocks.Types.SimpleController.PI or
          controllerTypeBoi == Modelica.Blocks.Types.SimpleController.PID,
          tab="Control", group="Boiler"));
  parameter Modelica.Units.SI.Time TdBoi(min=0)=10
    "Time constant of Derivative block" annotation (Dialog(enable=
          controllerTypeBoi == Modelica.Blocks.Types.SimpleController.PD or
          controllerTypeBoi == Modelica.Blocks.Types.SimpleController.PID,
          tab="Control", group="Boiler"));
  parameter Real wpBoi(min=0) = 1 "Set-point weight for Proportional block (0..1)"
    annotation (Dialog(tab="Control", group="Boiler"));
  parameter Real wdBoi(min=0) = 0 "Set-point weight for Derivative block (0..1)"
    annotation(Dialog(enable=
      controllerTypeBoi==.Modelica.Blocks.Types.SimpleController.PD or
      controllerTypeBoi==.Modelica.Blocks.Types.SimpleController.PID,
      tab="Control", group="Boiler"));
  parameter Real NiBoi(min=100*Modelica.Constants.eps) = 0.9
    "Ni*Ti is time constant of anti-windup compensation"
     annotation(Dialog(enable=
       controllerTypeBoi==.Modelica.Blocks.Types.SimpleController.PI or
       controllerTypeBoi==.Modelica.Blocks.Types.SimpleController.PID,
       tab="Control", group="Boiler"));
  parameter Real NdBoi(min=100*Modelica.Constants.eps) = 10
    "The higher Nd, the more ideal the derivative block"
    annotation(Dialog(enable=
      controllerTypeBoi==.Modelica.Blocks.Types.SimpleController.PD or
      controllerTypeBoi==.Modelica.Blocks.Types.SimpleController.PID,
      tab="Control", group="Boiler"));

  // Feedwater pump controller
  parameter Modelica.Blocks.Types.SimpleController controllerTypePum=
    Modelica.Blocks.Types.SimpleController.PI "Type of controller"
    annotation (Dialog(tab="Control", group="Pump"));
  parameter Real kPum(min=0) = 5 "Gain of controller"
    annotation (Dialog(tab="Control", group="Pump"));
  parameter Modelica.Units.SI.Time TiPum(min=Modelica.Constants.small)=120
    "Time constant of Integrator block"
    annotation (Dialog(enable=
      controllerTypePum == Modelica.Blocks.Types.SimpleController.PI or
      controllerTypePum == Modelica.Blocks.Types.SimpleController.PID,
      tab="Control", group="Pump"));
  parameter Modelica.Units.SI.Time TdPum(min=0)=0.1
    "Time constant of Derivative block"
    annotation (Dialog(enable=
      controllerTypePum == Modelica.Blocks.Types.SimpleController.PD or
      controllerTypePum == Modelica.Blocks.Types.SimpleController.PID,
      tab="Control", group="Pump"));
  parameter Real wpPum(min=0) = 1 "Set-point weight for Proportional block (0..1)"
    annotation (Dialog(tab="Control", group="Pump"));
  parameter Real wdPum(min=0) = 0 "Set-point weight for Derivative block (0..1)"
    annotation(Dialog(enable=
      controllerTypePum==.Modelica.Blocks.Types.SimpleController.PD or
      controllerTypePum==.Modelica.Blocks.Types.SimpleController.PID,
      tab="Control", group="Pump"));
  parameter Real NiPum(min=100*Modelica.Constants.eps) = 0.9
    "Ni*Ti is time constant of anti-windup compensation"
     annotation(Dialog(enable=
       controllerTypePum==.Modelica.Blocks.Types.SimpleController.PI or
       controllerTypePum==.Modelica.Blocks.Types.SimpleController.PID,
       tab="Control", group="Pump"));
  parameter Real NdPum(min=100*Modelica.Constants.eps) = 10
    "The higher Nd, the more ideal the derivative block";

 // Initial conditions
  parameter Modelica.Units.SI.Volume VTanFW_start=1
    "Setpoint for liquid water volume in the boiler"
    annotation(Dialog(tab = "Initialization"));
  parameter Modelica.Media.Interfaces.Types.AbsolutePressure pBoi_start=pSteSet
    "Start value of boiler pressure"
    annotation(Dialog(tab = "Initialization"));
  parameter Real yPum_start=0.7 "Initial value of output"
    annotation(Dialog(tab="Initialization"));

  // Setpoints
  parameter Modelica.Units.SI.Pressure pSteSet=600000
    "Steam pressure setpoint";

  parameter Modelica.Units.SI.AbsolutePressure pTanFW=101325
    "Pressure of feedwater tank";
  parameter Modelica.Units.SI.Volume VBoiWatSet=VBoi/2
    "Setpoint for liquid water volume in the boiler";

      // System sizing
  parameter Modelica.Units.SI.Volume VBoi=3
    "Total drum volume of steam boiler";

     // Dynamics
  parameter Modelica.Fluid.Types.Dynamics energyDynamics=
    Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
    "Type of energy balance: dynamic (3 initialization options) or steady state"
    annotation (Evaluate=true,Dialog(tab="Dynamics",group="Equations"));
  parameter Modelica.Fluid.Types.Dynamics massDynamics=energyDynamics
    "Type of mass balance: dynamic (3 initialization options) or steady state"
    annotation (Evaluate=true,Dialog(tab="Dynamics",group="Equations"));

  Buildings.Experimental.DHC.Loads.Steam.BaseClasses.ControlVolumeCondensation vol(
    redeclare final package MediumSte = MediumSte,
    redeclare final package MediumWat = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    final p_start=pSat,
    final m_flow_nominal=m_flow_nominal,
    V=1)
    "Volume"
    annotation (Placement(transformation(extent={{20,50},{40,30}})));
  Buildings.Fluid.FixedResistances.PressureDrop res(
    redeclare final package Medium = MediumWat,
    final m_flow_nominal=m_flow_nominal,
    final dp_nominal = dpPip)
    "Resistance in district network"
    annotation (Placement(transformation(extent={{0,-50},{-20,-30}})));
  Buildings.Experimental.DHC.Loads.Steam.BaseClasses.SteamTrap steTra(
    redeclare final package Medium = MediumWat,
    final m_flow_nominal=m_flow_nominal)
    "Steam trap"
    annotation (Placement(transformation(extent={{60,30},{80,50}})));
  Modelica.Blocks.Sources.Sine inp(
    amplitude=-0.5,
    f=1/86400,
    phase=3.1415926535898,
    offset=0.5)
    "Input signal"
    annotation (Placement(transformation(extent={{160,-20},{140,0}})));
  Buildings.Fluid.Sensors.MassFlowRate senMasFlo(
    redeclare final package Medium = MediumWat)
    "Mass flow rate sensor"
    annotation (Placement(transformation(extent={{80,-50},{60,-30}})));
  Buildings.Controls.Continuous.LimPID conPumCNR(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=1,
    Ti=15)
    "Controller"
    annotation (Placement(transformation(extent={{80,-20},{60,0}})));
  Modelica.Blocks.Math.Gain m_flow(final k=m_flow_nominal)
    "Gain to calculate m_flow"
    annotation (Placement(transformation(extent={{120,-20},{100,0}})));
  Buildings.Fluid.Movers.SpeedControlled_y pumCNR(
    redeclare final package Medium = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    final per=perPumCNR,
    y_start=1)
    "Condensate return pump"
    annotation (Placement(transformation(extent={{40,-50},{20,-30}})));
  Components.BoilerPolynomialFurnaceHeatBalance boi(
    m1_flow_nominal=1,
    m2_flow_nominal=m_flow_nominal,
    p_start=600000,
    redeclare package MediumAir = MediumAir,
    redeclare package MediumWat = MediumWat,
    redeclare package MediumSte = MediumSte,
    Q_flow_nominal=Q_flow_nominal,
    T_nominal=T_nominal,
    fue=Buildings.Fluid.Data.Fuels.NaturalGasHigherHeatingValue(),
    FA_ratio=1.15)
    annotation (Placement(transformation(extent={{-140,-6},{-120,14}})));
  Buildings.Controls.Continuous.LimPID conPum(
    final controllerType=Modelica.Blocks.Types.SimpleController.PI,
    final k=kPum,
    final Ti=TiPum,
    final Td=TdPum,
    final wp=wpPum,
    final wd=wdPum,
    final Ni=NiPum,
    final Nd=NdPum)
    "Pump control"
    annotation (Placement(transformation(extent={{-280,82},{-260,62}})));
  Modelica.Blocks.Sources.Constant uni(final k=1) "Unitary"
    annotation (Placement(transformation(extent={{-320,62},{-300,82}})));
  Fluid.Storage.ExpansionVessel           tanFW(
    redeclare package Medium = MediumWat,
    final V_start=VTanFW_start,
    final p_start=pTanFW)
    "Feedwater tank"
    annotation (Placement(transformation(extent={{-280,14},{-260,34}})));
  Fluid.Movers.SpeedControlled_y           pumFW(
    redeclare package Medium = MediumWat,
    final energyDynamics=energyDynamics,
    per=perPumFW,
    final y_start=yPum_start)
    "Feedwater pump"
    annotation (Placement(transformation(extent={{-240,-10},{-220,10}})));
  Fluid.FixedResistances.CheckValve           cheVal(
    redeclare package Medium = MediumWat,
    final m_flow_nominal=m_flow_nominal,
    dpValve_nominal=6000)
    "Check valve"
    annotation (Placement(transformation(extent={{-180,-10},{-160,10}})));
  Buildings.Controls.Continuous.LimPID conBoi(
    final controllerType=Modelica.Blocks.Types.SimpleController.PI,
    final k=kBoi,
    final Ti=TiBoi,
    final Td=TdBoi,
    final wp=wpBoi,
    final wd=wdBoi,
    final Ni=NiBoi,
    final Nd=NdBoi)
    "Boiler control"
    annotation (Placement(transformation(extent={{-200,80},{-180,60}})));
  Modelica.Blocks.Math.Gain PNor(final k=1/pSteSet)
    "Normalized pressure setpoint"
    annotation (Placement(transformation(extent={{-80,110},{-100,130}})));
  Fluid.Sensors.Pressure           senPreSte(redeclare final package
      Medium =
        MediumSte)
    "Steam pressure sensor"
    annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=-90,
        origin={-60,60})));
  Fluid.Sources.Boundary_pT exh(
    redeclare package Medium = MediumAir,
    T(displayUnit="K"),
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-60,12},{-80,32}})));
  Fluid.Sources.Boundary_pT pro(
    redeclare package Medium = MediumAir,
    use_T_in=false,
    T=573.15,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-202,10},{-182,30}})));
equation
  connect(vol.port_b, steTra.port_a)
    annotation (Line(points={{40,40},{60,40}}, color={0,127,255}));
  connect(steTra.port_b, senMasFlo.port_a) annotation (Line(points={{80,40},{90,
          40},{90,-40},{80,-40}}, color={0,127,255}));
  connect(senMasFlo.port_b, pumCNR.port_a)
    annotation (Line(points={{60,-40},{40,-40}},color={0,127,255}));
  connect(pumCNR.port_b, res.port_a)
    annotation (Line(points={{20,-40},{0,-40}},    color={0,127,255}));
  connect(senMasFlo.m_flow, conPumCNR.u_m)
    annotation (Line(points={{70,-29},{70,-22}}, color={0,0,127}));
  connect(conPumCNR.y, pumCNR.y)
    annotation (Line(points={{59,-10},{30,-10},{30,-28}},   color={0,0,127}));
  connect(m_flow.y, conPumCNR.u_s)
    annotation (Line(points={{99,-10},{82,-10}}, color={0,0,127}));
  connect(inp.y, m_flow.u)
    annotation (Line(points={{139,-10},{122,-10}},
                                                 color={0,0,127}));
  connect(pumFW.port_b, cheVal.port_a)
    annotation (Line(points={{-220,0},{-180,0}}, color={0,127,255}));
  connect(cheVal.port_b, boi.port_a2)
    annotation (Line(points={{-160,0},{-140,0}}, color={0,127,255}));
  connect(boi.port_b2, vol.port_a) annotation (Line(points={{-120,0},{0,0},{0,40},
          {20,40}}, color={0,127,255}));
  connect(senPreSte.port, boi.port_b2) annotation (Line(points={{-50,60},{-44,60},
          {-44,0},{-120,0}}, color={0,127,255}));
  connect(conBoi.y, boi.y) annotation (Line(points={{-179,70},{-150,70},{-150,12},
          {-142,12}}, color={0,0,127}));
  connect(senPreSte.p, PNor.u)
    annotation (Line(points={{-60,71},{-60,120},{-78,120}}, color={0,0,127}));
  connect(PNor.y, conBoi.u_m) annotation (Line(points={{-101,120},{-190,120},{-190,
          82}}, color={0,0,127}));
  connect(uni.y, conPum.u_s)
    annotation (Line(points={{-299,72},{-282,72}}, color={0,0,127}));
  connect(tanFW.port_a, pumFW.port_a)
    annotation (Line(points={{-270,14},{-270,0},{-240,0}}, color={0,127,255}));
  connect(tanFW.port_a, res.port_b) annotation (Line(points={{-270,14},{-270,-40},
          {-20,-40}}, color={0,127,255}));
  connect(conPum.y, pumFW.y)
    annotation (Line(points={{-259,72},{-230,72},{-230,12}}, color={0,0,127}));
  connect(boi.VLiq, conPum.u_m) annotation (Line(points={{-119,-4},{-102,-4},{-102,
          90},{-270,90},{-270,84}}, color={0,0,127}));
  connect(pro.ports[1], boi.port_a1) annotation (Line(points={{-182,20},{-152,20},
          {-152,8},{-140,8}}, color={0,127,255}));
  connect(boi.port_b1, exh.ports[1]) annotation (Line(points={{-120,8},{-86,8},{
          -86,22},{-80,22}}, color={0,127,255}));
  connect(uni.y, conBoi.u_s) annotation (Line(points={{-299,72},{-292,72},{-292,
          98},{-202,98},{-202,70}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-300,-100},
            {100,100}})),                                        Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-300,-100},{100,100}})),
    experiment(StopTime=86400, Tolerance=1e-6),
      __Dymola_Commands(file=
    "modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Plants/Steam/Examples/SingleBoiler.mos"
    "Simulate and plot"),
    Documentation(info="<html>
<p>This model validates the steam plant implemented in
<a href=\"modelica://Buildings.Experimental.DHC.Plants.Steam.SingleBoiler\">
Buildings.Experimental.DHC.Plants.Steam.SingleBoiler</a>.
</p>
</html>", revisions="<html>
<ul>
<li>
March 3, 2022 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"));
end SingleBoiler;
