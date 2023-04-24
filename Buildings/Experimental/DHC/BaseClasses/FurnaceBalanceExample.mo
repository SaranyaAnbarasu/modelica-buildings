within Buildings.Experimental.DHC.BaseClasses;
model FurnaceBalanceExample

  Buildings.GEDCalibration.CUBoulder.Components.FurnaceHeatBalance
    furnaceHeatBalance(
    T_exh_nominal(displayUnit="K") = 300.817,
    QFue=55533000,
    FA_ratio=1.05)
    annotation (Placement(transformation(extent={{-20,0},{0,20}})));
  Modelica.Blocks.Sources.Ramp           hAir(
    height=300750,
    duration=3000,
    offset=26130) "Exhaust temperature setpoint"
    annotation (Placement(transformation(extent={{-60,-4},{-40,16}})));
  Modelica.Blocks.Sources.RealExpression qLos(y=50000)
    "Exhaust temperature setpoint"
    annotation (Placement(transformation(extent={{-60,16},{-40,36}})));
equation
  connect(hAir.y, furnaceHeatBalance.ha)
    annotation (Line(points={{-39,6},{-22,6}}, color={0,0,127}));
  connect(qLos.y, furnaceHeatBalance.qLos) annotation (Line(points={{-39,26},{
          -28,26},{-28,16},{-22,16}}, color={0,0,127}));
  annotation ();
end FurnaceBalanceExample;
