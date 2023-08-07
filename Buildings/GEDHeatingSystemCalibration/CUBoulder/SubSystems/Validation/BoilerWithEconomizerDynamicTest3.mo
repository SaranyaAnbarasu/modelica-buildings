within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model BoilerWithEconomizerDynamicTest3
  extends Modelica.Icons.Example;

  //data file
  parameter String data = ("modelica://Buildings/Resources/Data/GEDCalibration/BoilerClosedLoopCase23328000.mos");

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
  parameter Modelica.Units.SI.AbsolutePressure p_nominal=900000
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

   // pumps
  parameter Buildings.Fluid.Movers.Data.Generic perPumFW(pressure(V_flow=
          m_flow_nominal*1000*{0.4,0.6,0.8,1.0}, dp=(p_nominal - 101325)*{1.34,1.27,
          1.17,1.0})) "Performance data for feedwater pump";
  parameter Buildings.Fluid.Movers.Data.Generic perPumCNR(pressure(V_flow=
          m_flow_nominal*1000*{0,1,2}, dp=dpPip*{2,1,0}))
    "Performance data for condensate return pumps";

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
  parameter Modelica.Media.Interfaces.Types.AbsolutePressure pBoi_start=p_nominal
    "Start value of boiler pressure"
    annotation(Dialog(tab = "Initialization"));
  parameter Real yPum_start=0.7 "Initial value of output"
    annotation(Dialog(tab="Initialization"));

//setpoints
 parameter Modelica.Units.SI.Volume VBoiWatSet=V/2
    "Setpoint for liquid water volume in the boiler";

    //buildings

  parameter Modelica.Units.SI.AbsolutePressure pSatBui=300000
    "Saturation pressure";
  parameter Modelica.Units.SI.AbsolutePressure pTanFW=101325
    "Pressure of feedwater tank";

      //resistances
  parameter Modelica.Units.SI.PressureDifference dpPip=6000
    "Pressure drop in the condensate return pipe";

  Fluid.Movers.FlowControlled_m_flow fwPum(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    per=perPumFW,
    addPowerToMedium=false,
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
    annotation (Placement(transformation(extent={{-416,-38},{-396,-18}})));
  Fluid.Sources.Boundary_pT pro(
    redeclare package Medium = MediumAir,
    use_T_in=false,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-180,-150},{-160,-130}})));
  HeatTransfer.Sources.FixedTemperature TAmb(T(displayUnit="K") = 294)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{160,-120},{140,-100}})));
  Fluid.Sources.Boundary_pT exh(
    redeclare package Medium = MediumAir,
    p=100000,
    T(displayUnit="K") = 560,
    nPorts=1) "Source" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={100,102})));

  Modelica.Blocks.Sources.Constant ecoValSig(k=1)
    annotation (Placement(transformation(extent={{-140,90},{-120,110}})));
  Fluid.Sensors.TemperatureTwoPort senTem(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={100,70})));
  Fluid.Sensors.TemperatureTwoPort senTem1(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={100,-100})));
  Fluid.Sensors.TemperatureTwoPort senTem2(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={106,-148})));
  Components.BoilerPolynomialFurnaceHeatBalance boi(
    m1_flow_nominal=m1_flow_nominal,
    m2_flow_nominal=m2_flow_nominal,
    show_T=true,
    energyDynamics=energyDynamics,
    massDynamics=massDynamics,
    p_start=p_nominal,
    redeclare package MediumAir = MediumAir,
    redeclare package MediumWat = MediumWat,
    redeclare package MediumSte = MediumSte,
    mDry=mDry,
    m_flow_nominal=10,
    dp_nominal=400000,
    Q_flow_nominal=Q_flow_nominal,
    fue=Buildings.Fluid.Data.Fuels.NaturalGasHigherHeatingValue(),
    UA=0.01*Q_flow_nominal/100,
    V=V,
    V_com=V_com,
    FA_ratio=FA_ratio,
    T_exh_nominal(displayUnit="K") = T_exh_nominal)
    annotation (Placement(transformation(extent={{12,-154},{32,-134}})));
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
    eps_nominal=0.6,
    r_nominal=1)   annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={94,-10})));

  Fluid.Actuators.Valves.TwoWayLinear val(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{40,10},{60,30}})));
  Fluid.Actuators.Valves.TwoWayLinear val1(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
    annotation (Placement(transformation(extent={{60,-50},{40,-30}})));
  Fluid.Actuators.Valves.TwoWayLinear val2(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=dpValve_nominal)
                           annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-26,-10})));
  Modelica.Blocks.Sources.Constant valSig(k=0)
    annotation (Placement(transformation(extent={{-140,-50},{-120,-30}})));

  Modelica.Blocks.Sources.CombiTimeTable T_ma(
    tableOnFile=true,
    table=fill(
        0.0,
        0,
        2),
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={16},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    "temperature of mixed air"
    annotation (Placement(transformation(extent={{-140,-110},{-120,-90}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    mFloFW1(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={2},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-468,-114},{-448,-94}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                    outputs(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns=2:15,
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-470,-72},{-450,-52}})));
  Modelica.Blocks.Sources.CombiTimeTable T_FW(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={3},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-428,4},{-408,24}})));
  Fluid.FixedResistances.PressureDrop res1(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m1_flow_nominal,
    dp_nominal=10000)
    annotation (Placement(transformation(extent={{140,-158},{160,-138}})));
  Fluid.Actuators.Dampers.PressureIndependent maDam(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m1_flow_nominal,
    dpDamper_nominal=1000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-100,-140})));
  Fluid.Actuators.Dampers.PressureIndependent staDam(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m1_flow_nominal,
    dpDamper_nominal=1000) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={100,32})));
  Fluid.Actuators.Dampers.PressureIndependent fgrDam(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m1_flow_nominal,
    dpDamper_nominal=1000) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={50,-80})));
  Fluid.Sensors.TemperatureTwoPort senTem3(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-12,-140})));
  Modelica.Blocks.Sources.CombiTimeTable yMa(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={13},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-140,-82},{-120,-62}})));
  Modelica.Blocks.Sources.CombiTimeTable ySta(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={15},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{160,22},{140,42}})));
  Modelica.Blocks.Sources.CombiTimeTable yFgr(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={14},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{160,-60},{140,-40}})));
  Fluid.Sensors.TemperatureTwoPort senTem4(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={10,-80})));
  Fluid.Sensors.TemperatureTwoPort senTem5(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K")) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-60,-140})));
  Fluid.HeatExchangers.Heater_T hea(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m1_flow_nominal,
    dp_nominal=1000)
    annotation (Placement(transformation(extent={{-140,-150},{-120,-130}})));
  Modelica.Blocks.Math.Gain m_fw_nom(k=6.30)
    annotation (Placement(transformation(extent={{-96,54},{-84,66}})));
  Controls.FeedwaterSingleElementControl feedwaterSingleElementControl(
      controllerType=Modelica.Blocks.Types.SimpleController.PI, VBoiWatSet=
        VBoiWatSet)
    annotation (Placement(transformation(extent={{-140,50},{-120,70}})));
  Experimental.DHC.Loads.Steam.BaseClasses.ValveSelfActing prv(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m1_flow_nominal,
    pb_nominal=pSatBui,
    dp_start(displayUnit="bar") = 600000)
    annotation (Placement(transformation(extent={{180,-158},{200,-138}})));
  Experimental.DHC.Loads.Steam.BaseClasses.ControlVolumeCondensation vol(
    redeclare final package MediumSte = MediumSte,
    redeclare final package MediumWat = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    final p_start=pSatBui,
    final m_flow_nominal=m1_flow_nominal,
    V=10) "Volume"
    annotation (Placement(transformation(extent={{222,-138},{242,-158}})));
  Experimental.DHC.Loads.Steam.BaseClasses.SteamTrap steTra(redeclare final
      package Medium = MediumWat, final m_flow_nominal=m1_flow_nominal)
    "Steam trap" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={262,-178})));
  Fluid.Sensors.MassFlowRate senMasFlo(redeclare final package Medium =
        MediumWat) "Mass flow rate sensor"
    annotation (Placement(transformation(extent={{242,-190},{222,-210}})));
  Fluid.Movers.SpeedControlled_y pumCNR(
    redeclare final package Medium = MediumWat,
    energyDynamics=energyDynamics,
    final per=perPumCNR,
    y_start=1) "Condensate return pump"
    annotation (Placement(transformation(extent={{202,-210},{182,-190}})));
  Buildings.Controls.Continuous.LimPID conPumCNR(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=1,
    Ti=15) "Controller"
    annotation (Placement(transformation(extent={{316,-178},{296,-158}})));
  Modelica.Blocks.Sources.CombiTimeTable mFloFW2(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns={6},
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{364,-178},{344,-158}})));
  Fluid.FixedResistances.PressureDrop res(
    redeclare final package Medium = MediumWat,
    final m_flow_nominal=m1_flow_nominal,
    final dp_nominal=dpPip) "Resistance in district network"
    annotation (Placement(transformation(extent={{-88,-210},{-108,-190}})));
  Controls.Firing firing(pBoiSet=p_nominal)
    annotation (Placement(transformation(extent={{-70,-80},{-50,-60}})));
  Fluid.Sensors.Pressure           senPreSte(redeclare final package Medium =
        MediumSte)
    "Steam pressure sensor"
    annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=0,
        origin={58,-168})));
  Fluid.FixedResistances.CheckValve           cheVal(
    redeclare package Medium = MediumWat,
    final m_flow_nominal=m1_flow_nominal,
    dpValve_nominal=6000)
    "Check valve"
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-26,-64})));
  Fluid.Sources.Boundary_pT sou(redeclare package Medium = MediumWat, nPorts=2)
              "Source" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-206,18})));
  Fluid.FixedResistances.PressureDrop res2(
    redeclare final package Medium = MediumWat,
    final m_flow_nominal=m1_flow_nominal,
    final dp_nominal=dpPip) "Resistance in district network"
    annotation (Placement(transformation(extent={{-140,10},{-120,30}})));
equation

  connect(eco.port_a2, val.port_b) annotation (Line(
      points={{88,0},{88,20},{60,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_a, eco.port_b2) annotation (Line(
      points={{60,-40},{88,-40},{88,-20}},
      color={0,127,255},
      thickness=0.5));
  connect(val2.port_a, fwPum.port_b) annotation (Line(
      points={{-26,0},{-26,20},{-50,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val.port_a, fwPum.port_b) annotation (Line(
      points={{40,20},{-50,20}},
      color={0,127,255},
      thickness=0.5));
  connect(TAmb.port, boi.heatPort)
    annotation (Line(points={{140,-110},{22,-110},{22,-134}},
                                                          color={191,0,0}));
  connect(boi.port_b1, senTem1.port_a) annotation (Line(
      points={{32,-140},{100,-140},{100,-110}},
      color={0,0,0},
      thickness=0.5));
  connect(boi.port_b2, senTem2.port_a) annotation (Line(
      points={{32,-148},{96,-148}},
      color={238,46,47},
      thickness=0.5));
  connect(exh.ports[1], senTem.port_b) annotation (Line(
      points={{100,92},{100,80}},
      color={0,0,0},
      thickness=0.5));
  connect(ecoValSig.y, val.y) annotation (Line(
      points={{-119,100},{50,100},{50,32}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig.y, val1.y) annotation (Line(
      points={{-119,100},{70,100},{70,-6},{50,-6},{50,-28}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(valSig.y, val2.y) annotation (Line(
      points={{-119,-40},{-44,-40},{-44,-10},{-38,-10}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senTem2.port_b, res1.port_a) annotation (Line(
      points={{116,-148},{140,-148}},
      color={238,46,47},
      thickness=0.5));
  connect(senTem3.port_a, boi.port_a1)
    annotation (Line(points={{-2,-140},{12,-140}}, color={141,137,137},
      thickness=0.5));
  connect(senTem.port_a, staDam.port_b)
    annotation (Line(points={{100,60},{100,42}}, color={0,0,0},
      thickness=0.5));
  connect(staDam.port_a, eco.port_b1)
    annotation (Line(points={{100,22},{100,0}}, color={0,0,0},
      thickness=0.5));
  connect(eco.port_a1, senTem1.port_b)
    annotation (Line(points={{100,-20},{100,-90}}, color={0,0,0},
      thickness=0.5));
  connect(yMa.y[1], maDam.y) annotation (Line(points={{-119,-72},{-100,-72},{-100,
          -128}},                     color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ySta.y[1], staDam.y)
    annotation (Line(points={{139,32},{112,32}}, color={0,0,127},
      pattern=LinePattern.Dash));
  connect(yFgr.y[1], fgrDam.y)
    annotation (Line(points={{139,-50},{50,-50},{50,-68}}, color={0,0,127},
      pattern=LinePattern.Dash));
  connect(fgrDam.port_a, senTem1.port_b) annotation (Line(points={{60,-80},{100,
          -80},{100,-90}}, color={0,0,0},
      thickness=0.5));
  connect(fgrDam.port_b, senTem4.port_a)
    annotation (Line(points={{40,-80},{20,-80}}, color={0,0,0},
      thickness=0.5));
  connect(senTem4.port_b, senTem3.port_b) annotation (Line(points={{0,-80},{-34,
          -80},{-34,-140},{-22,-140}}, color={0,0,0},
      thickness=0.5));
  connect(maDam.port_b, senTem5.port_b)
    annotation (Line(points={{-90,-140},{-70,-140}}, color={141,137,137},
      thickness=0.5));
  connect(senTem5.port_a, senTem3.port_b)
    annotation (Line(points={{-50,-140},{-22,-140}}, color={141,137,137},
      thickness=0.5));
  connect(pro.ports[1], hea.port_a)
    annotation (Line(points={{-160,-140},{-140,-140}}, color={141,137,137},
      thickness=0.5));
  connect(hea.port_b, maDam.port_a)
    annotation (Line(points={{-120,-140},{-110,-140}}, color={141,137,137},
      thickness=0.5));
  connect(T_ma.y[1], hea.TSet) annotation (Line(points={{-119,-100},{-114,-100},
          {-114,-124},{-148,-124},{-148,-132},{-142,-132}},
                                    color={0,0,127},
      pattern=LinePattern.Dash));
  connect(m_fw_nom.y, fwPum.m_flow_in)
    annotation (Line(points={{-83.4,60},{-60,60},{-60,32}},
                                                          color={0,0,127},
      pattern=LinePattern.Dash));
  connect(feedwaterSingleElementControl.y, m_fw_nom.u) annotation (Line(
      points={{-119,60},{-97.2,60}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(boi.VLiq, feedwaterSingleElementControl.VWatBoi) annotation (Line(
      points={{33,-152},{40,-152},{40,-100},{-10,-100},{-10,80},{-148,80},{-148,
          60},{-142,60}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senMasFlo.m_flow,conPumCNR. u_m) annotation (Line(
      points={{232,-211},{232,-216},{306,-216},{306,-180}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(conPumCNR.y,pumCNR. y) annotation (Line(
      points={{295,-168},{192,-168},{192,-188}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(conPumCNR.u_s,mFloFW2. y[1])
    annotation (Line(points={{318,-168},{343,-168}}, color={0,0,127},
      pattern=LinePattern.Dash));
  connect(pumCNR.port_b,res. port_a)
    annotation (Line(points={{182,-200},{-88,-200}},color={0,127,255}));
  connect(res1.port_b, prv.port_a)
    annotation (Line(points={{160,-148},{180,-148}}, color={0,127,255}));
  connect(prv.port_b, vol.port_a)
    annotation (Line(points={{200,-148},{222,-148}}, color={0,127,255}));
  connect(vol.port_b, steTra.port_a) annotation (Line(points={{242,-148},{262,-148},
          {262,-168}}, color={0,127,255}));
  connect(steTra.port_b, senMasFlo.port_a) annotation (Line(points={{262,-188},{
          262,-200},{242,-200}}, color={0,127,255}));
  connect(senMasFlo.port_b, pumCNR.port_a)
    annotation (Line(points={{222,-200},{202,-200}}, color={0,127,255}));
  connect(firing.y, boi.y) annotation (Line(
      points={{-49,-70},{-42,-70},{-42,-124},{4,-124},{4,-136},{10,-136}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senPreSte.port, boi.port_b2) annotation (Line(
      points={{58,-158},{58,-148},{32,-148}},
      color={238,46,47},
      thickness=0.5));
  connect(senPreSte.p, firing.PBoiLin) annotation (Line(
      points={{47,-168},{-80,-168},{-80,-70},{-72,-70}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(val2.port_b, cheVal.port_a) annotation (Line(
      points={{-26,-20},{-26,-54}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_b, cheVal.port_a) annotation (Line(
      points={{40,-40},{-26,-40},{-26,-54}},
      color={0,127,255},
      thickness=0.5));
  connect(cheVal.port_b, boi.port_a2) annotation (Line(
      points={{-26,-74},{-26,-148},{12,-148}},
      color={0,127,255},
      thickness=0.5));
  connect(res.port_b, sou.ports[1]) annotation (Line(
      points={{-108,-200},{-190,-200},{-190,19},{-196,19}},
      color={0,127,255},
      thickness=0.5));
  connect(sou.ports[2], res2.port_a) annotation (Line(
      points={{-196,17},{-194,17},{-194,20},{-140,20}},
      color={0,127,255},
      thickness=0.5));
  connect(res2.port_b, fwPum.port_a) annotation (Line(
      points={{-120,20},{-70,20}},
      color={0,127,255},
      thickness=0.5));
  annotation (
    __Dymola_Commands(file="modelica://Buildings/Resources/Data/GEDCalibration/Scripts/BoilerWithEconomizerClosedLoopImproved.mos"
        "Simulate and plot"),
    experiment(
      StartTime=23328000,
      StopTime=24000019.2,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,
            100}})));
end BoilerWithEconomizerDynamicTest3;
