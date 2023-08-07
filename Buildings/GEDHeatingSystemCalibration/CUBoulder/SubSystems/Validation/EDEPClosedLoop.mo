within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model EDEPClosedLoop
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

    //economizer
      parameter Real r_nominal=2/3
    "Ratio between air-side and water-side convective heat transfer (hA-value) at nominal condition";
  parameter Modelica.Units.SI.Temperature T_a2_nominal(displayUnit="K") = 250
    "Nominal temperature at port a2";
  parameter Modelica.Units.SI.Temperature T_a1_nominal(displayUnit="K") = 577
    "Nominal temperature at port a1";

      parameter Modelica.Units.SI.MassFlowRate m1_flow_nominal=3.43
    "Nominal mass flow rate";

  parameter Modelica.Units.SI.MassFlowRate m2_flow_nominal=6.36
    "Nominal mass flow rate";

  parameter Modelica.Units.SI.PressureDifference dp1_nominal = 20000
    "Pressure drop at m_flow_nominal";

   parameter Modelica.Units.SI.PressureDifference dp2_nominal = 39989.6
    "Pressure drop at m_flow_nominal";

  Fluid.Movers.FlowControlled_m_flow fwPum(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m_flow_nominal,
    addPowerToMedium=false,
    nominalValuesDefineDefaultPressureCurve=true,
    dp_nominal=dp_nominal) "Feed water pump"
    annotation (Placement(transformation(extent={{-140,10},{-120,30}})));
  Modelica.Blocks.Sources.Trapezoid
                               mFloFW(
    amplitude=5.9,
    rising=300,
    width=300,
    falling=300,
    period=900)
    annotation (Placement(transformation(extent={{-240,106},{-220,126}})));
  Modelica.Blocks.Sources.SawTooth maDamSig(
    amplitude=80,
    period=2000,
    offset=273)
    annotation (Placement(transformation(extent={{-240,-22},{-220,-2}})));
  Fluid.Sources.Boundary_pT pro(
    redeclare package Medium = MediumAir,
    use_T_in=true,
    T=573.15,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-160,-102},{-140,-82}})));
  HeatTransfer.Sources.FixedTemperature           TAmb(T(displayUnit="K") = 284)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-50,70})));
  Fluid.Sources.Boundary_pT exh(
    redeclare package Medium = MediumAir,
    p=100000,
    T(displayUnit="K") = 560,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={66,70})));

  Modelica.Blocks.Sources.Constant ecoValSig(k=1)
    annotation (Placement(transformation(extent={{-240,40},{-220,60}})));
  Fluid.Sensors.TemperatureTwoPort
                                 senTem(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K"))
    annotation (Placement(transformation(extent={{-10,10},{10,-10}},
        rotation=90,
        origin={66,38})));
  Fluid.Sensors.TemperatureTwoPort
                                 senTem2(
    redeclare package Medium = MediumSte,
    m_flow_nominal=m_flow_nominal,
    tau=30,
    T_start(displayUnit="K"))
    annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=180,
        origin={90,-100})));
  Modelica.Blocks.Sources.Constant
                               Tfw_in(k=427)
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-230,-88})));
  Components.BoilerPolynomialFurnaceHeatBalance
                                     boi(
    m1_flow_nominal=10,
    m2_flow_nominal=m_flow_nominal,
    show_T=true,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    p_start=1000000,
    redeclare package MediumAir = MediumAir,
    redeclare package MediumWat = MediumWat,
    redeclare package MediumSte = MediumSte,
    mDry=1.5E-3*Q_flow_nominal,
    m_flow_nominal=10,
    dp_nominal=400000,
    Q_flow_nominal=Q_flow_nominal,
    fue(
      h=46402971,
      d=800,
      mCO2=2.2),
    UA=0.01*Q_flow_nominal/100,
    V=50,
    V_com=20,
    FA_ratio=1.15,
    T_exh_nominal(displayUnit="K") = 420)
    annotation (Placement(transformation(extent={{12,-106},{32,-86}})));
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
    use_Q_flow_nominal=true,
    Q_flow_nominal=1324000,
    T_a1_nominal(displayUnit="K") = 677,
    T_a2_nominal(displayUnit="K") = 250,
    r_nominal=0.1)
    annotation (Placement(transformation(extent={{-10,10},{10,-10}},
        rotation=90,
        origin={60,4})));

  Fluid.Actuators.Valves.TwoWayLinear val(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m_flow_nominal,
    dpValve_nominal=50000)
    annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
  Fluid.Actuators.Valves.TwoWayLinear val1(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m_flow_nominal,
    dpValve_nominal=50000)
    annotation (Placement(transformation(extent={{-40,-26},{-60,-6}})));
  Fluid.Actuators.Valves.TwoWayLinear val2(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m_flow_nominal,
    dpValve_nominal=50000) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-92,2})));
  Modelica.Blocks.Sources.Constant ecoDamSig(k=1)
    annotation (Placement(transformation(extent={{140,-18},{120,2}})));
  Modelica.Blocks.Sources.Constant ecoValSig1(k=1)
    annotation (Placement(transformation(extent={{-240,70},{-220,90}})));
  Fluid.Storage.ExpansionVessel           tanFW(
    redeclare package Medium = MediumWat,
    final V_start=VTanFW_start,
    final p_start=pTanFW)
    "Feedwater tank"
    annotation (Placement(transformation(extent={{-190,38},{-170,58}})));
  Experimental.DHC.Loads.Steam.BaseClasses.ControlVolumeCondensation           vol(
    redeclare final package MediumSte = MediumSte,
    redeclare final package MediumWat = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    final p_start=pSat,
    final m_flow_nominal=m_flow_nominal,
    V=1)
    "Volume"
    annotation (Placement(transformation(extent={{120,-90},{140,-110}})));
  Fluid.FixedResistances.PressureDrop           res(
    redeclare final package Medium = MediumWat,
    final m_flow_nominal=m_flow_nominal,
    final dp_nominal=dpPip)
    "Resistance in district network"
    annotation (Placement(transformation(extent={{20,-190},{0,-170}})));
  Experimental.DHC.Loads.Steam.BaseClasses.SteamTrap           steTra(
      redeclare final package Medium = MediumWat, final m_flow_nominal=
        m_flow_nominal)
    "Steam trap"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={160,-142})));
  Modelica.Blocks.Sources.Sine inp(
    amplitude=-0.5,
    f=1/86400,
    phase=3.1415926535898,
    offset=0.5)
    "Input signal"
    annotation (Placement(transformation(extent={{280,-130},{260,-110}})));
  Fluid.Sensors.MassFlowRate           senMasFlo(redeclare final package
      Medium = MediumWat)
    "Mass flow rate sensor"
    annotation (Placement(transformation(extent={{140,-170},{120,-190}})));
  Buildings.Controls.Continuous.LimPID conPumCNR(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=1,
    Ti=15)
    "Controller"
    annotation (Placement(transformation(extent={{200,-130},{180,-110}})));
  Modelica.Blocks.Math.Gain m_flow(final k=m_flow_nominal)
    "Gain to calculate m_flow"
    annotation (Placement(transformation(extent={{236,-130},{216,-110}})));
  Fluid.Movers.SpeedControlled_y           pumCNR(
    redeclare final package Medium = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    final per=perPumCNR,
    y_start=1)
    "Condensate return pump"
    annotation (Placement(transformation(extent={{100,-190},{80,-170}})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd(redeclare package
      Medium = MediumAir) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={66,-30})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd1(redeclare package
              Medium = MediumAir) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={32,-60})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd2(redeclare package
              Medium = MediumAir) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-92})));
  Fluid.HeatExchangers.Heater_T hea annotation (Placement(transformation(
          extent={{-120,-102},{-100,-82}})));
  Modelica.Blocks.Sources.Constant fgrDamSig(k=1)
    annotation (Placement(transformation(extent={{140,-50},{120,-30}})));
  Modelica.Blocks.Sources.Constant
                               Tfw_in1(k=427)
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-230,-50})));
equation

  connect(eco.port_a2, val.port_b) annotation (Line(
      points={{54,14},{54,20},{-40,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_a, eco.port_b2) annotation (Line(
      points={{-40,-16},{54,-16},{54,-6}},
      color={244,125,35},
      thickness=0.5));
  connect(val2.port_a, fwPum.port_b) annotation (Line(
      points={{-92,12},{-92,20},{-120,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val2.port_b, boi.port_a2) annotation (Line(points={{-92,-8},{
          -92,-100},{12,-100}},
                      color={0,127,255}));
  connect(val.port_a, fwPum.port_b) annotation (Line(
      points={{-60,20},{-120,20}},
      color={0,127,255},
      thickness=0.5));
  connect(TAmb.port, boi.heatPort)
    annotation (Line(points={{-40,70},{0,70},{0,-80},{22,-80},{22,-86}},
                                                          color={191,0,0}));
  connect(boi.port_b2, senTem2.port_a) annotation (Line(
      points={{32,-100},{80,-100}},
      color={238,46,47},
      thickness=0.5));
  connect(exh.ports[1], senTem.port_b) annotation (Line(
      points={{66,60},{66,48}},
      color={0,0,0},
      thickness=0.5));
  connect(senTem.port_a, eco.port_b1) annotation (Line(
      points={{66,28},{66,14}},
      color={0,0,0},
      thickness=0.5));
  connect(tanFW.port_a, fwPum.port_a) annotation (Line(
      points={{-180,38},{-180,20},{-140,20}},
      color={0,127,255},
      thickness=0.5));
  connect(vol.port_b,steTra. port_a)
    annotation (Line(points={{140,-100},{160,-100},{160,-132}},
                                               color={238,46,47},
      thickness=0.5));
  connect(steTra.port_b,senMasFlo. port_a) annotation (Line(points={{160,
          -152},{160,-180},{140,-180}},
                                  color={0,127,255},
      thickness=0.5));
  connect(senMasFlo.port_b,pumCNR. port_a)
    annotation (Line(points={{120,-180},{100,-180}},
                                                color={0,127,255},
      thickness=0.5));
  connect(pumCNR.port_b,res. port_a)
    annotation (Line(points={{80,-180},{20,-180}}, color={0,127,255},
      thickness=0.5));
  connect(senMasFlo.m_flow,conPumCNR. u_m)
    annotation (Line(points={{130,-191},{130,-196},{190,-196},{190,-132}},
                                                 color={0,0,127},
      pattern=LinePattern.Dash));
  connect(conPumCNR.y,pumCNR. y)
    annotation (Line(points={{179,-120},{90,-120},{90,-168}},
                                                            color={0,0,127},
      pattern=LinePattern.Dash));
  connect(m_flow.y,conPumCNR. u_s)
    annotation (Line(points={{215,-120},{202,-120}},
                                                 color={0,0,127},
      pattern=LinePattern.Dash));
  connect(inp.y,m_flow. u)
    annotation (Line(points={{259,-120},{238,-120}},
                                                 color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senTem2.port_b, vol.port_a) annotation (Line(
      points={{100,-100},{120,-100}},
      color={238,46,47},
      thickness=0.5));
  connect(val2.port_b, val1.port_b) annotation (Line(
      points={{-92,-8},{-92,-16},{-60,-16}},
      color={0,127,255},
      thickness=0.5));
  connect(boi.port_b1, damPreInd.port_a) annotation (Line(
      points={{32,-92},{66,-92},{66,-40}},
      color={0,0,0},
      thickness=0.5));
  connect(eco.port_a1, damPreInd.port_b) annotation (Line(
      points={{66,-6},{66,-20}},
      color={0,0,0},
      thickness=0.5));
  connect(damPreInd1.port_b, damPreInd.port_a) annotation (Line(
      points={{42,-60},{66,-60},{66,-40}},
      color={0,0,0},
      thickness=0.5));
  connect(damPreInd1.port_a, boi.port_a1) annotation (Line(
      points={{22,-60},{-40,-60},{-40,-92},{12,-92}},
      color={0,0,0},
      thickness=0.5));
  connect(pro.ports[1], hea.port_a) annotation (Line(
      points={{-140,-92},{-120,-92}},
      color={144,142,142},
      thickness=0.5));
  connect(hea.port_b, damPreInd2.port_a) annotation (Line(
      points={{-100,-92},{-80,-92}},
      color={144,142,142},
      thickness=0.5));
  connect(tanFW.port_a, res.port_b) annotation (Line(
      points={{-180,38},{-180,-180},{0,-180}},
      color={0,127,255},
      thickness=0.5));
  connect(ecoDamSig.y, damPreInd.y) annotation (Line(
      points={{119,-8},{60,-8},{60,-30},{54,-30}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(fgrDamSig.y, damPreInd1.y) annotation (Line(
      points={{119,-40},{32,-40},{32,-48}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig.y, val2.y) annotation (Line(
      points={{-219,50},{-200,50},{-200,2},{-104,2}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig1.y, val.y) annotation (Line(
      points={{-219,80},{-66,80},{-66,40},{-50,40},{-50,32}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig1.y, val1.y) annotation (Line(
      points={{-219,80},{-66,80},{-66,4},{-50,4},{-50,-4}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(mFloFW.y, fwPum.m_flow_in) annotation (Line(
      points={{-219,116},{-130,116},{-130,32}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(Tfw_in.y, pro.T_in) annotation (Line(
      points={{-219,-88},{-162,-88}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(Tfw_in1.y, hea.TSet) annotation (Line(
      points={{-219,-50},{-128,-50},{-128,-84},{-122,-84}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(maDamSig.y, damPreInd2.y) annotation (Line(
      points={{-219,-12},{-118,-12},{-118,-74},{-70,-74},{-70,-80}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(val1.port_b, boi.port_a2) annotation (Line(
      points={{-60,-16},{-92,-16},{-92,-100},{12,-100}},
      color={244,125,35},
      thickness=0.5));
  connect(damPreInd2.port_b, boi.port_a1) annotation (Line(
      points={{-60,-92},{12,-92}},
      color={144,142,142},
      thickness=0.5));
  annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Components/HeatBalanceBoiler.mos"
        "Simulate and plot"),
    experiment(Tolerance=1e-6, StopTime=3600),Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
            {160,100}})),                                        Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})));
end EDEPClosedLoop;
