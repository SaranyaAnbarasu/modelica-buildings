within Buildings.Experimental.DHC.Plants.Steam.BaseClasses;
model Combustion "Empirical model for combustion process"
  extends Buildings.Fluid.Interfaces.PartialTwoPortInterface(allowFlowReversal=false);

 Modelica.Blocks.Interfaces.RealInput y(min=0, max=1) "Part load ratio"
    annotation (Placement(transformation(extent={{-140,50},{-100,90}}),
        iconTransformation(extent={{-120,70},{-100,90}})));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-92,4},{92,-4}},
          lineColor={244,125,35},
          lineThickness=0.5,
          fillColor={244,125,35},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-80,82},{80,-60}},
          lineColor={244,125,35},
          lineThickness=0.5,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-38,42},{38,-36}},
          lineColor={238,46,47},
          lineThickness=0.5,
          fillColor={238,46,47},
          fillPattern=FillPattern.Sphere)}),                     Diagram(
        coordinateSystem(preserveAspectRatio=false)));








end Combustion;
