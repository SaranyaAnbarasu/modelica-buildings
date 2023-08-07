within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model HeatBalanceBoilerWithEconomizerClosedLoopDynamic
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
    annotation (Placement(transformation(extent={{-100,10},{-80,30}})));
  Modelica.Blocks.Sources.CombiTimeTable
                               mFloFW
    annotation (Placement(transformation(extent={{-200,90},{-180,110}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                   maDamSig(
    offset=273)
    annotation (Placement(transformation(extent={{-200,-40},{-180,-20}})));
  Fluid.Sources.Boundary_pT pro(
    redeclare package Medium = MediumAir,
    use_T_in=true,
    T=573.15,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{-138,-102},{-118,-82}})));
  HeatTransfer.Sources.FixedTemperature           TAmb(T(displayUnit="K") = 284)
    "Ambient temperature in boiler room"
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-10,70})));
  Fluid.Sources.Boundary_pT exh(
    redeclare package Medium = MediumAir,
    p=100000,
    T(displayUnit="K") = 560,
    nPorts=1) "Source"
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={66,70})));

  Modelica.Blocks.Sources.Constant ecoValSig(k=1)
    annotation (Placement(transformation(extent={{-200,-8},{-180,12}})));
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
  Modelica.Blocks.Sources.CombiTimeTable
                               Tfw_in
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-190,-88})));
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
    annotation (Placement(transformation(extent={{22,-106},{42,-86}})));
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
    annotation (Placement(transformation(extent={{-20,10},{0,30}})));
  Fluid.Actuators.Valves.TwoWayLinear val1(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m_flow_nominal,
    dpValve_nominal=50000)
    annotation (Placement(transformation(extent={{0,-30},{-20,-10}})));
  Fluid.Actuators.Valves.TwoWayLinear val2(
    redeclare package Medium = MediumWat,
    m_flow_nominal=m_flow_nominal,
    dpValve_nominal=50000) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-60,2})));
  Modelica.Blocks.Sources.CombiTimeTable
                                   ecoDamSig
    annotation (Placement(transformation(extent={{140,-18},{120,2}})));
  Modelica.Blocks.Sources.Constant ecoValSig1(k=1)
    annotation (Placement(transformation(extent={{-200,60},{-180,80}})));
  Fluid.Storage.ExpansionVessel           tanFW(
    redeclare package Medium = MediumWat,
    final V_start=VTanFW_start,
    final p_start=pTanFW)
    "Feedwater tank"
    annotation (Placement(transformation(extent={{-170,38},{-150,58}})));
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
    annotation (Placement(transformation(extent={{20,-170},{0,-150}})));
  Experimental.DHC.Loads.Steam.BaseClasses.SteamTrap           steTra(
      redeclare final package Medium = MediumWat, final m_flow_nominal=
        m_flow_nominal)
    "Steam trap"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={160,-130})));
  Modelica.Blocks.Sources.CombiTimeTable
                               inp(
    offset=0.5)
    "Input signal"
    annotation (Placement(transformation(extent={{242,-120},{222,-100}})));
  Fluid.Sensors.MassFlowRate           senMasFlo(redeclare final package
      Medium = MediumWat)
    "Mass flow rate sensor"
    annotation (Placement(transformation(extent={{140,-150},{120,-170}})));
  Buildings.Controls.Continuous.LimPID conPumCNR(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=1,
    Ti=15)
    "Controller"
    annotation (Placement(transformation(extent={{200,-120},{180,-100}})));
  Fluid.Movers.SpeedControlled_y           pumCNR(
    redeclare final package Medium = MediumWat,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    final per=perPumCNR,
    y_start=1)
    "Condensate return pump"
    annotation (Placement(transformation(extent={{100,-170},{80,-150}})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd(redeclare package
      Medium = MediumAir) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={66,-30})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd1(redeclare package
              Medium = MediumAir) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={30,-60})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd2(redeclare package
              Medium = MediumAir) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-30,-92})));
  Fluid.HeatExchangers.Heater_T hea annotation (Placement(transformation(
          extent={{-98,-102},{-78,-82}})));
  Modelica.Blocks.Sources.CombiTimeTable
                                   fgrDamSig
    annotation (Placement(transformation(extent={{140,-50},{120,-30}})));
  Modelica.Blocks.Sources.CombiTimeTable
                               Tfw_in1
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-190,-60})));
equation

  connect(eco.port_a2, val.port_b) annotation (Line(
      points={{54,14},{54,20},{0,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val1.port_a, eco.port_b2) annotation (Line(
      points={{0,-20},{54,-20},{54,-6}},
      color={244,125,35},
      thickness=0.5));
  connect(val2.port_a, fwPum.port_b) annotation (Line(
      points={{-60,12},{-60,20},{-80,20}},
      color={0,127,255},
      thickness=0.5));
  connect(val2.port_b, boi.port_a2) annotation (Line(points={{-60,-8},{
          -60,-100},{22,-100}},
                      color={0,127,255}));
  connect(val.port_a, fwPum.port_b) annotation (Line(
      points={{-20,20},{-80,20}},
      color={0,127,255},
      thickness=0.5));
  connect(TAmb.port, boi.heatPort)
    annotation (Line(points={{0,70},{10,70},{10,-80},{32,-80},{32,-86}},
                                                          color={191,0,0}));
  connect(boi.port_b2, senTem2.port_a) annotation (Line(
      points={{42,-100},{80,-100}},
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
      points={{-160,38},{-160,20},{-100,20}},
      color={0,127,255},
      thickness=0.5));
  connect(vol.port_b,steTra. port_a)
    annotation (Line(points={{140,-100},{160,-100},{160,-120}},
                                               color={238,46,47},
      thickness=0.5));
  connect(steTra.port_b,senMasFlo. port_a) annotation (Line(points={{160,
          -140},{160,-160},{140,-160}},
                                  color={0,127,255},
      thickness=0.5));
  connect(senMasFlo.port_b,pumCNR. port_a)
    annotation (Line(points={{120,-160},{100,-160}},
                                                color={0,127,255},
      thickness=0.5));
  connect(pumCNR.port_b,res. port_a)
    annotation (Line(points={{80,-160},{20,-160}}, color={0,127,255},
      thickness=0.5));
  connect(senMasFlo.m_flow,conPumCNR. u_m)
    annotation (Line(points={{130,-171},{130,-178},{190,-178},{190,-122}},
                                                 color={0,0,127},
      pattern=LinePattern.Dash));
  connect(conPumCNR.y,pumCNR. y)
    annotation (Line(points={{179,-110},{90,-110},{90,-148}},
                                                            color={0,0,127},
      pattern=LinePattern.Dash));
  connect(senTem2.port_b, vol.port_a) annotation (Line(
      points={{100,-100},{120,-100}},
      color={238,46,47},
      thickness=0.5));
  connect(val2.port_b, val1.port_b) annotation (Line(
      points={{-60,-8},{-60,-20},{-20,-20}},
      color={0,127,255},
      thickness=0.5));
  connect(boi.port_b1, damPreInd.port_a) annotation (Line(
      points={{42,-92},{66,-92},{66,-40}},
      color={0,0,0},
      thickness=0.5));
  connect(eco.port_a1, damPreInd.port_b) annotation (Line(
      points={{66,-6},{66,-20}},
      color={0,0,0},
      thickness=0.5));
  connect(damPreInd1.port_b, damPreInd.port_a) annotation (Line(
      points={{40,-60},{66,-60},{66,-40}},
      color={0,0,0},
      thickness=0.5));
  connect(damPreInd1.port_a, boi.port_a1) annotation (Line(
      points={{20,-60},{0,-60},{0,-92},{22,-92}},
      color={0,0,0},
      thickness=0.5));
  connect(pro.ports[1], hea.port_a) annotation (Line(
      points={{-118,-92},{-98,-92}},
      color={144,142,142},
      thickness=0.5));
  connect(hea.port_b, damPreInd2.port_a) annotation (Line(
      points={{-78,-92},{-40,-92}},
      color={144,142,142},
      thickness=0.5));
  connect(tanFW.port_a, res.port_b) annotation (Line(
      points={{-160,38},{-160,-160},{0,-160}},
      color={0,127,255},
      thickness=0.5));
  connect(ecoDamSig.y, damPreInd.y) annotation (Line(
      points={{119,-8},{60,-8},{60,-30},{54,-30}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(fgrDamSig.y, damPreInd1.y) annotation (Line(
      points={{119,-40},{30,-40},{30,-48}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig.y, val2.y) annotation (Line(
      points={{-179,2},{-72,2}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig1.y, val.y) annotation (Line(
      points={{-179,70},{-40,70},{-40,40},{-10,40},{-10,32}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(ecoValSig1.y, val1.y) annotation (Line(
      points={{-179,70},{-40,70},{-40,4},{-10,4},{-10,-8}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(mFloFW.y, fwPum.m_flow_in) annotation (Line(
      points={{-179,100},{-90,100},{-90,32}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(Tfw_in.y, pro.T_in) annotation (Line(
      points={{-179,-88},{-140,-88}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(Tfw_in1.y, hea.TSet) annotation (Line(
      points={{-179,-60},{-106,-60},{-106,-84},{-100,-84}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(maDamSig.y, damPreInd2.y) annotation (Line(
      points={{-179,-30},{-80,-30},{-80,-74},{-30,-74},{-30,-80}},
      color={0,0,127},
      pattern=LinePattern.Dash));
  connect(val1.port_b, boi.port_a2) annotation (Line(
      points={{-20,-20},{-60,-20},{-60,-100},{22,-100}},
      color={244,125,35},
      thickness=0.5));
  connect(damPreInd2.port_b, boi.port_a1) annotation (Line(
      points={{-20,-92},{22,-92}},
      color={144,142,142},
      thickness=0.5));
  connect(inp.y[1], conPumCNR.u_s)
    annotation (Line(points={{221,-110},{202,-110}}, color={0,0,127}));
  annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Components/HeatBalanceBoiler.mos"
        "Simulate and plot"),
    experiment(Tolerance=1e-6, StopTime=3600),Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
            {160,100}})),                                        Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})));
end HeatBalanceBoilerWithEconomizerClosedLoopDynamic;
