within Buildings.GEDHeatingSystemCalibration.CUBoulder.Controls;
block Firing

  Modelica.Blocks.Sources.Constant uni(final k=1) "Unitary"
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  Buildings.Controls.Continuous.LimPID conPum(
    final controllerType=controllerType,
    final k=k,
    final Ti=Ti,
    final Td=Td,
    final wp=wp,
    final wd=wd,
    final Ni=Ni,
    final Nd=Nd)
    "Pump control"
    annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=180,
        origin={10,0})));
  Modelica.Blocks.Math.Gain PNor(final k=1/pBoiSet)
    "Normalized pressure setpoint" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-30,-40})));
  Modelica.Blocks.Interfaces.RealInput PBoiLin "Boiler headder line pressure"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}}),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
          extent={{100,-10},{120,10}}), iconTransformation(extent={{100,-10},{120,
            10}})));
  parameter Real k=10
                    "Gain of controller";
  parameter Modelica.Units.SI.Time Ti=120 "Time constant of Integrator block";
  parameter Real wp=1 "Set-point weight for Proportional block (0..1)";
  parameter Real Ni=0.9 "Ni*Ti is time constant of anti-windup compensation";
  parameter Modelica.Units.SI.Time Td=10 "Time constant of Derivative block";
  parameter Modelica.Blocks.Types.SimpleController controllerType=Modelica.Blocks.Types.SimpleController.PI
    "Type of controller";
  parameter Real wd=0 "Set-point weight for Derivative block (0..1)";
  parameter Real Nd=10 "The higher Nd, the more ideal the derivative block";

  //setpoints
  parameter Modelica.Units.SI.AbsolutePressure pBoiSet = 900000 "Boiler headder pressure setpoint";

equation
  connect(conPum.u_s, uni.y)
    annotation (Line(points={{-2,0},{-39,0}}, color={0,0,127}));
  connect(conPum.u_m, PNor.y)
    annotation (Line(points={{10,-12},{10,-40},{-19,-40}}, color={0,0,127}));
  connect(PNor.u,PBoiLin)
    annotation (Line(points={{-42,-40},{-82,-40},{-82,0},{-120,0}},
                                                    color={0,0,127}));
  connect(conPum.y, y) annotation (Line(points={{21,-8.88178e-16},{52,-8.88178e-16},
          {52,0},{110,0}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
  Diagram(coordinateSystem(preserveAspectRatio=false)));
end Firing;
