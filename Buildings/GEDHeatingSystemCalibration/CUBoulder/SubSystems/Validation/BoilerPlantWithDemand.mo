within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model BoilerPlantWithDemand
  extends Modelica.Icons.Example;

  //data file
  parameter String data=("modelica://Buildings/Resources/Data/GEDCalibration/BoilerWithEconCase23328000.mos");

  // Medium declarations
  package MediumWat =
      Buildings.Media.Specialized.Water.TemperatureDependentDensity (p_default=101325,
        T_default=273.15 + 100) "Water medium - port_a (inlet)";
  package MediumSte = Buildings.Media.Steam (
      p_default=300000,
      T_default=273.15 + 200,
      h_default=2700000) "Steam medium - port_b (oulet)";
  package MediumAir = Buildings.Media.CombustionAir "Combustion air medium";

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
  parameter Modelica.Media.Interfaces.Types.AbsolutePressure p_start=900000
    "Start value of pressure";

  //economizer
  parameter Real r_nominal=2/3
    "Ratio between air-side and water-side convective heat transfer (hA-value) at nominal condition";
  parameter Modelica.Units.SI.Power Q_eco_flow_nominal=1324000
    "Nominal power of economizer";

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

  //resistances
  parameter Modelica.Units.SI.PressureDifference dpPip=6000
    "Pressure drop in the condensate return pipe";

  // pumps
  parameter Buildings.Fluid.Movers.Data.Generic perPumFW(pressure(V_flow=
          m_flow_nominal*1000*{0.4,0.6,0.8,1.0}, dp=(pSteSet - 101325)*{1.34,1.27,
          1.17,1.0})) "Performance data for feedwater pump";
  parameter Buildings.Fluid.Movers.Data.Generic perPumCNR(pressure(V_flow=
          m_flow_nominal*1000*{0,1,2}, dp=dpPip*{2,1,0}))
    "Performance data for condensate return pumps";

  //buildings

  parameter Modelica.Units.SI.AbsolutePressure pSatBui=300000
    "Saturation pressure";
  parameter Modelica.Units.SI.AbsolutePressure pTanFW=101325
    "Pressure of feedwater tank";

  parameter Modelica.Units.SI.Volume VTanFW_start=10
    "Setpoint for liquid water volume in the boiler"
    annotation (Dialog(tab="Initialization"));

    parameter Modelica.Units.SI.PressureDifference dpDamper_nominal = 100
    "Pressure drop of fully open damper at nominal mass flow rate";

  // Setpoints
  parameter Modelica.Units.SI.Pressure pSteSet=900000
    "Steam pressure setpoint";

      parameter Modelica.Units.SI.Volume VBoiWatSet=V/2
    "Setpoint for liquid water volume in the boiler";

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

  parameter Real yPum_start=0.7 "Initial value of output"
    annotation(Dialog(tab="Initialization"));

  // Dynamics
  parameter Modelica.Fluid.Types.Dynamics energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
    "Type of energy balance: dynamic (3 initialization options) or steady state"
    annotation (Evaluate=true, Dialog(tab="Dynamics", group="Equations"));
  parameter Modelica.Fluid.Types.Dynamics massDynamics=energyDynamics
    "Type of mass balance: dynamic (3 initialization options) or steady state"
    annotation (Evaluate=true, Dialog(tab="Dynamics", group="Equations"));

  Fluid.Movers.SpeedControlled_y     fwPum(
    redeclare package Medium = MediumWat,
    per=perPumFW,
    addPowerToMedium=false)
                           "Feed water pump"
    annotation (Placement(transformation(extent={{-70,10},{-50,30}})));
  Modelica.Blocks.Sources.CombiTimeTable plr(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={4},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-408,-190},{-388,-170}})));
  Modelica.Blocks.Sources.Constant Tair_in(k=294.15)
    annotation (Placement(transformation(extent={{-392,4},{-372,24}})));
  Fluid.Sources.Boundary_pT pro(
    redeclare package Medium = MediumAir,
    use_T_in=false,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-150,-142},{-130,-122}})));
  HeatTransfer.Sources.FixedTemperature TAmb(T(displayUnit="K") = 294)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{160,-110},{140,-90}})));
  Fluid.Sources.Boundary_pT exh(
    redeclare package Medium = MediumAir,
    T(displayUnit="K"),
    nPorts=1) "Source" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={66,70})));

  Modelica.Blocks.Sources.Constant ecoValSig(k=1)
    annotation (Placement(transformation(extent={{-220,70},{-200,90}})));
  Fluid.Sensors.TemperatureTwoPort senTem(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={66,38})));
  Fluid.Sensors.TemperatureTwoPort senTem2(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={72,-140})));
  Modelica.Blocks.Sources.Constant Tfw_in(k=377) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-386,-90})));
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
    annotation (Placement(transformation(extent={{10,-146},{30,-126}})));

  Fluid.Actuators.Valves.TwoWayLinear val(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{10,10},{30,30}})));
  Fluid.Actuators.Valves.TwoWayLinear val1(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{30,-38},{10,-18}})));
  Fluid.Actuators.Valves.TwoWayLinear val2(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-20,0})));
  Modelica.Blocks.Sources.Constant valSig(k=1)
    annotation (Placement(transformation(extent={{-220,-10},{-200,10}})));

  Modelica.Blocks.Sources.Sine inp(
    amplitude=-0.5,
    f=1/86400,
    phase=3.1415926535898,
    offset=0.5) "Input signal"
    annotation (Placement(transformation(extent={{-396,62},{-376,82}})));
  Modelica.Blocks.Sources.CombiTimeTable T_ma(
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
    "Temperature of Mixed air"
    annotation (Placement(transformation(extent={{-332,-128},{-312,-108}})));
  Modelica.Blocks.Sources.CombiTimeTable mFloFW1(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={2},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-400,-158},{-380,-138}})));
  Modelica.Blocks.Sources.CombiTimeTable outputs(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns=2:12,
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-360,-50},{-340,-30}})));
  Modelica.Blocks.Sources.CombiTimeTable T_FW(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={5},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-360,-20},{-340,0}})));
  Fluid.FixedResistances.PressureDrop res1(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m1_flow_nominal,
    dp_nominal=10000)
    annotation (Placement(transformation(extent={{98,-150},{118,-130}})));
  Modelica.Blocks.Sources.Constant ecoDamSig(k=0.5)
    annotation (Placement(transformation(extent={{-346,36},{-366,56}})));
  Modelica.Blocks.Sources.Constant fgrDamSig(k=0.5)
    annotation (Placement(transformation(extent={{-346,68},{-366,88}})));
  Fluid.HeatExchangers.Heater_T hea(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m1_flow_nominal,
    dp_nominal=1000)
    annotation (Placement(transformation(extent={{-120,-142},{-100,-122}})));
  Fluid.FixedResistances.PressureDrop res(
    redeclare final package Medium = MediumWat,
    final m_flow_nominal=m1_flow_nominal,
    final dp_nominal=dpPip) "Resistance in district network"
    annotation (Placement(transformation(extent={{32,-210},{12,-190}})));
  Fluid.Movers.SpeedControlled_y pumCNR(
    redeclare final package Medium = MediumWat,
    energyDynamics=energyDynamics,
    final per=perPumCNR,
    y_start=1) "Condensate return pump"
    annotation (Placement(transformation(extent={{160,-210},{140,-190}})));
  Fluid.Sensors.MassFlowRate senMasFlo(redeclare final package Medium =
        MediumWat) "Mass flow rate sensor"
    annotation (Placement(transformation(extent={{200,-190},{180,-210}})));
  Experimental.DHC.Loads.Steam.BaseClasses.SteamTrap steTra(redeclare final
      package       Medium =
                       MediumWat, final m_flow_nominal=m1_flow_nominal)
    "Steam trap" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={220,-170})));
  Buildings.Controls.Continuous.LimPID conPumCNR(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=1,
    Ti=15) "Controller"
    annotation (Placement(transformation(extent={{274,-170},{254,-150}})));
  Experimental.DHC.Loads.Steam.BaseClasses.ControlVolumeCondensation vol(
    redeclare final package MediumSte = MediumSte,
    redeclare final package MediumWat = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    final p_start=pSatBui,
    final m_flow_nominal=m1_flow_nominal,
    V=10) "Volume"
    annotation (Placement(transformation(extent={{180,-130},{200,-150}})));
  Fluid.Storage.ExpansionVessel tanFW(
    redeclare package Medium = MediumWat,
    final V_start=VTanFW_start,
    final p_start=pTanFW) "Feedwater tank"
    annotation (Placement(transformation(extent={{-170,20},{-150,40}})));
  Modelica.Blocks.Sources.CombiTimeTable mFloFW2(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={2},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{322,-170},{302,-150}})));
  Modelica.Blocks.Sources.Constant fgrDamSig1(k=1)
    annotation (Placement(transformation(extent={{-398,-222},{-378,-202}})));

  Experimental.DHC.Loads.Steam.BaseClasses.ValveSelfActing prv(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m1_flow_nominal,
    pb_nominal=pSatBui,
    dp_start=p_nominal - pSatBui)
    annotation (Placement(transformation(extent={{138,-150},{158,-130}})));
  Fluid.Sensors.TemperatureTwoPort senTem1(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-24,-132})));
  Fluid.FixedResistances.CheckValve           cheVal(
    redeclare package Medium = MediumWat,
    final m_flow_nominal=m_flow_nominal,
    dpValve_nominal=6000)
    "Check valve"
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-20,-62})));
  Fluid.Sensors.Pressure           senPreSte(redeclare final package
      Medium =
        MediumSte)
    "Steam pressure sensor"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={56,-170})));
  Modelica.Blocks.Math.Gain PNor(final k=1/pSteSet)
    "Normalized pressure setpoint"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-70,-104})));
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
    annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
  Modelica.Blocks.Sources.Constant uni(final k=1) "Unitary"
    annotation (Placement(transformation(extent={{-126,-60},{-106,-40}})));
  Modelica.Blocks.Math.Gain PNor1(final k=1/V)
    "Normalized pressure setpoint"
    annotation (Placement(transformation(extent={{-14,104},{-34,124}})));
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
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-60,114})));
  Modelica.Blocks.Sources.Constant uni1(final k=1)
                                                  "Unitary"
    annotation (Placement(transformation(extent={{-124,118},{-104,138}})));
  Modelica.Blocks.Sources.Constant steHtr(k=400)
    annotation (Placement(transformation(extent={{-214,-114},{-194,-94}})));
equation

  connect(val2.port_a, fwPum.port_b) annotation (Line(
      points={{-20,10},{-20,20},{-50,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val.port_a, fwPum.port_b) annotation (Line(
      points={{10,20},{-50,20}},
      color={0,127,255},
      thickness=0.5));
  connect(TAmb.port, boi.heatPort) annotation (Line(points={{140,-100},{20,-100},
          {20,-126}}, color={191,0,0}));
  connect(boi.port_b2, senTem2.port_a) annotation (Line(
      points={{30,-140},{62,-140}},
      color={238,46,47},
      thickness=0.5));
  connect(exh.ports[1], senTem.port_b) annotation (Line(
      points={{66,60},{66,48}},
      color={0,0,0},
      thickness=0.5));
  connect(ecoValSig.y, val.y) annotation (Line(
      points={{-199,80},{0,80},{0,38},{20,38},{20,32}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig.y, val1.y) annotation (Line(
      points={{-199,80},{0,80},{0,-10},{20,-10},{20,-16}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(valSig.y, val2.y) annotation (Line(
      points={{-199,0},{-196,0},{-196,6.66134e-16},{-32,6.66134e-16}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senTem2.port_b, res1.port_a) annotation (Line(
      points={{82,-140},{98,-140}},
      color={238,46,47},
      thickness=0.5));
  connect(conPumCNR.y, pumCNR.y) annotation (Line(
      points={{253,-160},{150,-160},{150,-188}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senMasFlo.m_flow, conPumCNR.u_m) annotation (Line(
      points={{190,-211},{190,-216},{264,-216},{264,-172}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(fwPum.port_a, tanFW.port_a)
    annotation (Line(points={{-70,20},{-160,20}}, color={0,127,255}));
  connect(pro.ports[1], hea.port_a)
    annotation (Line(points={{-130,-132},{-120,-132}}, color={0,127,255}));
  connect(vol.port_b, steTra.port_a) annotation (Line(points={{200,-140},
          {220,-140},{220,-160}},
                       color={0,127,255}));
  connect(steTra.port_b, senMasFlo.port_a) annotation (Line(points={{220,
          -180},{220,-200},{200,-200}},
                                 color={0,127,255}));
  connect(senMasFlo.port_b, pumCNR.port_a)
    annotation (Line(points={{180,-200},{160,-200}}, color={0,127,255}));
  connect(pumCNR.port_b, res.port_a)
    annotation (Line(points={{140,-200},{32,-200}}, color={0,127,255}));
  connect(res.port_b, tanFW.port_a) annotation (Line(points={{12,-200},{-180,-200},
          {-180,20},{-160,20}}, color={0,127,255}));
  connect(conPumCNR.u_s, mFloFW2.y[1])
    annotation (Line(points={{276,-160},{301,-160}}, color={0,0,127}));
  connect(res1.port_b, prv.port_a)
    annotation (Line(points={{118,-140},{138,-140}}, color={0,127,255}));
  connect(prv.port_b, vol.port_a)
    annotation (Line(points={{158,-140},{180,-140}}, color={0,127,255}));
  connect(senTem1.port_a, boi.port_a1)
    annotation (Line(points={{-14,-132},{10,-132}}, color={0,127,255}));
  connect(hea.port_b, senTem1.port_b)
    annotation (Line(points={{-100,-132},{-34,-132}}, color={0,127,255}));
  connect(val2.port_b, cheVal.port_a)
    annotation (Line(points={{-20,-10},{-20,-52}}, color={0,127,255}));
  connect(val1.port_b, cheVal.port_a) annotation (Line(points={{10,-28},{-20,-28},
          {-20,-52}}, color={0,127,255}));
  connect(cheVal.port_b, boi.port_a2) annotation (Line(points={{-20,-72},{-14,-72},
          {-14,-152},{10,-152},{10,-140}}, color={0,127,255}));
  connect(senPreSte.port, senTem2.port_a) annotation (Line(points={{56,-160},{56,
          -140},{62,-140}}, color={0,127,255}));
  connect(PNor.y, conBoi.u_m)
    annotation (Line(points={{-70,-93},{-70,-62}}, color={0,0,127}));
  connect(PNor.u, senPreSte.p) annotation (Line(points={{-70,-116},{-70,-172},{45,
          -172},{45,-170}}, color={0,0,127}));
  connect(uni.y, conBoi.u_s)
    annotation (Line(points={{-105,-50},{-82,-50}}, color={0,0,127}));
  connect(conBoi.y, boi.y) annotation (Line(points={{-59,-50},{-46,-50},{-46,-128},
          {8,-128}}, color={0,0,127}));
  connect(boi.VLiq, PNor1.u) annotation (Line(points={{31,-144},{34,-144},{34,-138},
          {40,-138},{40,114},{-12,114}},           color={0,0,127}));
  connect(uni1.y, conPum.u_s)
    annotation (Line(points={{-103,128},{-103,126},{-60,126}},
                                                             color={0,0,127}));
  connect(PNor1.y, conPum.u_m) annotation (Line(points={{-35,114},{-48,114}},
                           color={0,0,127}));
  connect(conPum.y, fwPum.y)
    annotation (Line(points={{-60,103},{-60,32}}, color={0,0,127}));
  connect(steHtr.y, hea.TSet) annotation (Line(points={{-193,-104},{-122,-104},{
          -122,-124}}, color={0,0,127}));
  connect(val.port_b, val1.port_a) annotation (Line(points={{30,20},{48,
          20},{48,-28},{30,-28}}, color={0,127,255}));
  connect(senTem.port_a, boi.port_b1) annotation (Line(points={{66,28},{
          66,-132},{30,-132}}, color={0,127,255}));
  annotation (
    __Dymola_Commands(file="modelica://Buildings/Resources/Data/GEDCalibration/Scripts/BoilerWithEconomizerDynamic.mos"
        "Simulate and plot"),
    experiment(
      StopTime=3600,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,
            100}})));
end BoilerPlantWithDemand;
