within Buildings.GEDHeatingSystemCalibration;
package UTAustin

  package Components

    model BoilerPolynomialSuperheated
      "Example model for the steam boiler with a polynomial efficiency curve"
      extends Modelica.Icons.Example;

      // Medium declarations
      package MediumWat =
          Buildings.Media.Specialized.Water.TemperatureDependentDensity
        "Water medium - port_a (inlet)";
      package MediumSte = Buildings.Media.Steam
         "Steam medium - port_b (oulet)";

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

      Modelica.Blocks.Sources.TimeTable y(
        table=[0,0; 1200,1; 1200,0; 2000,0; 2000,1; 3600,1])
        "Load ratio"
        annotation (Placement(transformation(extent={{156,34},{176,54}})));
      Buildings.Fluid.Sources.Boundary_pT sin(
        redeclare package Medium = MediumSte,
        p(displayUnit="bar") = 300000,
        T=423.15,
        nPorts=1)
        "Sink"
        annotation (Placement(transformation(extent={{92,-20},{72,0}})));
      Buildings.Fluid.Sources.Boundary_pT sou(
        redeclare package Medium = MediumWat,
        p=300000 + dp_nominal,
        T=303.15,
        nPorts=1)
        "Source"
        annotation (Placement(transformation(extent={{-138,-20},{-118,0}})));
      Buildings.HeatTransfer.Sources.FixedTemperature TAmb(T=288.15)
        "Ambient temperature in boiler room"
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Buildings.Experimental.DHC.Plants.Steam.BaseClasses.BoilerPolynomial boiDyn(
        p_start=1100000,
        redeclare package MediumSte = MediumSte,
        redeclare package MediumWat = MediumWat,
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        m_flow_nominal=m_flow_nominal,
        Q_flow_nominal=Q_flow_nominal,
        fue=Buildings.Fluid.Data.Fuels.NaturalGasLowerHeatingValue(),
        dp_nominal=dp_nominal)
        "Steam boiler with dynamic balance"
        annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
      Modelica.Blocks.Math.Gain PNor(final k=1/pSteSet)
        "Normalized pressure setpoint"
        annotation (Placement(transformation(extent={{-20,-86},{-40,-66}})));
      Fluid.Sensors.Pressure           senPreSte(redeclare final package
          Medium =
            MediumHea_b)
        "Steam pressure sensor"
        annotation (Placement(transformation(extent={{30,-66},{10,-86}})));
      Buildings.Controls.Continuous.LimPID conBoi(
        final controllerType=controllerTypeBoi,
        final k=kBoi,
        final Ti=TiBoi,
        final Td=TdBoi,
        final wp=wpBoi,
        final wd=wdBoi,
        final Ni=NiBoi,
        final Nd=NdBoi) "Boiler control"
        annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
      Modelica.Blocks.Math.Gain VNor(final k=1/VBoiWatSet)
        "Normalized volume setpoint"
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            rotation=-90,
            origin={40,30})));
      Buildings.Controls.Continuous.LimPID conPum(
        final controllerType=controllerTypePum,
        final k=kPum,
        final Ti=TiPum,
        final Td=TdPum,
        final wp=wpPum,
        final wd=wdPum,
        final Ni=NiPum,
        final Nd=NdPum) "Pump control"
        annotation (Placement(transformation(extent={{-148,86},{-128,66}})));
      Fluid.Movers.SpeedControlled_y           pumFW(
        final energyDynamics=energyDynamics,
        redeclare final package Medium = Medium,
        final per=per,
        final y_start=yPum_start)
        "Feedwater pump"
        annotation (Placement(transformation(extent={{-90,-20},{-70,0}})));
      Modelica.Blocks.Sources.Constant uni(final k=1) "Unitary"
        annotation (Placement(transformation(extent={{-188,66},{-168,86}})));
      Modelica.Blocks.Sources.Constant uni1(final k=1)
                                                      "Unitary"
        annotation (Placement(transformation(extent={{-120,-60},{-100,-40}})));
    equation
      connect(TAmb.port, boiDyn.heatPort)
        annotation (Line(points={{-20,50},{-10,50},{-10,-2.8}},
                                                              color={191,0,0}));
      connect(boiDyn.port_b, sin.ports[1])
        annotation (Line(points={{0,-10},{72,-10}},         color={0,127,255}));
      connect(boiDyn.port_a, pumFW.port_b)
        annotation (Line(points={{-20,-10},{-70,-10}}, color={0,127,255}));
      connect(pumFW.port_a, sou.ports[1])
        annotation (Line(points={{-90,-10},{-118,-10}}, color={0,127,255}));
      connect(boiDyn.port_b, senPreSte.port)
        annotation (Line(points={{0,-10},{20,-10},{20,-66}}, color={0,127,255}));
      connect(senPreSte.p, PNor.u)
        annotation (Line(points={{9,-76},{-18,-76}}, color={0,0,127}));
      connect(conPum.u_s, uni.y)
        annotation (Line(points={{-150,76},{-167,76}}, color={0,0,127}));
      connect(boiDyn.VLiq, VNor.u)
        annotation (Line(points={{1,-4},{40,-4},{40,18}}, color={0,0,127}));
      connect(VNor.y, conPum.u_m) annotation (Line(points={{40,41},{40,106},{-138,
              106},{-138,88}}, color={0,0,127}));
      connect(conPum.y, pumFW.y) annotation (Line(points={{-127,76},{-120,76},{-120,
              28},{-80,28},{-80,2}}, color={0,0,127}));
      connect(uni1.y, conBoi.u_s)
        annotation (Line(points={{-99,-50},{-82,-50}}, color={0,0,127}));
      connect(PNor.y, conBoi.u_m)
        annotation (Line(points={{-41,-76},{-70,-76},{-70,-62}}, color={0,0,127}));
      connect(conBoi.y, boiDyn.y) annotation (Line(points={{-59,-50},{-52,-50},{-52,
              -2},{-22,-2}}, color={0,0,127}));
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
</html>",     revisions="<html>
<ul>
<li>
July 23, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"),
        Diagram(coordinateSystem(extent={{-240,-100},{200,100}})),
        Icon(coordinateSystem(extent={{-240,-100},{200,100}})));
    end BoilerPolynomialSuperheated;

    model BoilerPolynomialSuperheatedWoControls
      "Example model for the steam boiler with a polynomial efficiency curve"
      extends Modelica.Icons.Example;

      // Medium declarations
      package MediumWat =
          Buildings.Media.Specialized.Water.TemperatureDependentDensity
        "Water medium - port_a (inlet)";
      package MediumSte = Buildings.Media.Steam
         "Steam medium - port_b (oulet)";

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

      Modelica.Blocks.Sources.TimeTable y(
        table=[0,0; 1200,1; 1200,0; 2000,0; 2000,1; 3600,1]) "Load ratio"
        annotation (Placement(transformation(extent={{-120,60},{-100,80}})));
      Buildings.Fluid.Sources.Boundary_pT sin(
        redeclare package Medium = MediumSte,
        p(displayUnit="bar") = 1100000,
        T=461.15,
        nPorts=1)
        "Sink"
        annotation (Placement(transformation(extent={{40,-20},{20,0}})));
      Buildings.Fluid.Sources.Boundary_pT sou(
        redeclare package Medium = MediumWat,
        p=300000 + dp_nominal,
        T=303.15,
        nPorts=1)
        "Source"
        annotation (Placement(transformation(extent={{-138,-20},{-118,0}})));
      Buildings.HeatTransfer.Sources.FixedTemperature TAmb(T=288.15)
        "Ambient temperature in boiler room"
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Buildings.Experimental.DHC.Plants.Steam.BaseClasses.BoilerPolynomial boiDyn(
        p_start=1100000,
        redeclare package MediumSte = MediumSte,
        redeclare package MediumWat = MediumWat,
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        m_flow_nominal=m_flow_nominal,
        Q_flow_nominal=Q_flow_nominal,
        fue=Buildings.Fluid.Data.Fuels.NaturalGasLowerHeatingValue(),
        dp_nominal=dp_nominal)
        "Steam boiler with dynamic balance"
        annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
      Fluid.Movers.FlowControlled_m_flow       pumFW(
        redeclare final package Medium = MediumWat,
        m_flow_nominal=0.3,
        nominalValuesDefineDefaultPressureCurve=true,
        dp_nominal(displayUnit="bar") = 1100000)
        "Feedwater pump"
        annotation (Placement(transformation(extent={{-100,-20},{-80,0}})));
      Modelica.Blocks.Sources.TimeTable FW_mflo(table=[0,0.002; 1200,0.002; 1200,0.002;
            2000,0.002; 2000,0.002; 3600,0.002]) "Feed water mass flow rate"
        annotation (Placement(transformation(extent={{-120,22},{-100,42}})));
    equation
      connect(TAmb.port, boiDyn.heatPort)
        annotation (Line(points={{-20,50},{-10,50},{-10,-2.8}},
                                                              color={191,0,0}));
      connect(boiDyn.port_b, sin.ports[1])
        annotation (Line(points={{0,-10},{20,-10}},         color={0,127,255}));
      connect(boiDyn.port_a, pumFW.port_b)
        annotation (Line(points={{-20,-10},{-80,-10}}, color={0,127,255}));
      connect(pumFW.port_a, sou.ports[1])
        annotation (Line(points={{-100,-10},{-118,-10}},color={0,127,255}));
      connect(y.y, boiDyn.y) annotation (Line(points={{-99,70},{-54,70},{-54,-2},{
              -22,-2}}, color={0,0,127}));
      connect(FW_mflo.y, pumFW.m_flow_in)
        annotation (Line(points={{-99,32},{-90,32},{-90,2}}, color={0,0,127}));
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
</html>",     revisions="<html>
<ul>
<li>
July 23, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"),
        Diagram(coordinateSystem(extent={{-160,-100},{60,100}})),
        Icon(coordinateSystem(extent={{-160,-100},{60,100}})));
    end BoilerPolynomialSuperheatedWoControls;
  end Components;

  package Controls
  end Controls;
end UTAustin;
