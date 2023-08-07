within Buildings.GEDHeatingSystemCalibration.CUBoulder.SubSystems.Validation;
model CombiTableDataCheck
    //data file
  parameter String data = ("modelica://Buildings/Resources/Data/GEDCalibration/BoilerWithEconCase21600000.mos");

  Modelica.Blocks.Sources.CombiTimeTable
                                    outputs(
    tableOnFile=true,
    tableName="table",
    fileName=ModelicaServices.ExternalReferences.loadResource(data),
    verboseRead=true,
    columns=2:9)
    annotation (Placement(transformation(extent={{-56,28},{-36,48}})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m1_flow_nominal,
    dpDamper_nominal=dpDamper_nominal) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-64,-60})));
  Fluid.Actuators.Dampers.PressureIndependent damPreInd2(
    redeclare package Medium = MediumAir,
    m_flow_nominal=m1_flow_nominal,
    dpDamper_nominal=dpDamper_nominal) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={-18,-30})));
  annotation (                                 experiment(
      StartTime=21600000,
      StopTime=31968000,
      __Dymola_Algorithm="Dassl"));
end CombiTableDataCheck;
