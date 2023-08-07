within Buildings.GEDHeatingSystemCalibration.CUBoulder;
package Components

  model BoilerPolynomialExhaust "Existing steam boiler model with additional exhaust side modeled; does include a constant efficiency; A equilibrium boiler with water phase change from liquid to vapor, discharging
  saturated steam vapor, with the efficiency curve described by a polynomial."
    extends Buildings.BaseClasses.BaseIconLow;
    extends
      Buildings.Experimental.DHC.BaseClasses.Steam.PartialFourPortInterfaceThreeMedium(
      redeclare final package Medium1 = MediumAir,
      redeclare final package Medium2 = MediumWat,
      redeclare final package Medium3 = MediumSte,
      final T_start=MediumSte.saturationTemperature(p_start));

    // Medium declarations
    replaceable package MediumAir =
        Buildings.Media.CombustionAir
      "Air medium on the combustion side";
    replaceable package MediumWat =
        Buildings.Media.Specialized.Water.TemperatureDependentDensity
      "Water medium - port_a2 (inlet)";
    replaceable package MediumSte = Buildings.Media.Steam
      "Steam medium - port_b2 (oulet)";
    // Initialization
    parameter Boolean fixed_p_start=false "Set to true if p_start is to be used as an explicit initial equation, 
    not an initial guess"   annotation (Dialog(tab="Initialization"));
    // Nominal conditions
    parameter Modelica.Units.SI.PressureDifference dp_nominal(displayUnit="Pa")
      "Pressure drop at nominal mass flow rate"
      annotation (Dialog(group="Nominal condition"));
    parameter Modelica.Units.SI.Power Q_flow_nominal "Nominal heating power";
    parameter Modelica.Units.SI.Temperature T_nominal=373.15 "Temperature used to compute nominal efficiency 
    (only used if efficiency curve depends on temperature)";

    // Efficiency, fuel, and boiler properties
    parameter Buildings.Fluid.Types.EfficiencyCurves effCur=Buildings.Fluid.Types.EfficiencyCurves.Constant
      "Curve used to compute the efficiency";
    parameter Real a[:]={0.9} "Coefficients for efficiency curve";
    parameter Buildings.Fluid.Data.Fuels.Generic fue "Fuel type"
      annotation (choicesAllMatching=true);
    parameter Modelica.Units.SI.ThermalConductance UA=0.05*Q_flow_nominal/30
      "Overall UA value";
    parameter Modelica.Units.SI.Volume V=1.5E-6*Q_flow_nominal
      "Total internal volume of boiler" annotation (Dialog(tab="Dynamics", enable=
           not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

    parameter Modelica.Units.SI.Volume V_com=m1_flow_nominal*tau1/rho1_nominal
      "Total internal volume of combustion" annotation (Dialog(tab="Dynamics",
          enable=not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

    parameter Modelica.Units.SI.Mass mDry=1.5E-3*Q_flow_nominal
      "Mass of boiler that will be lumped to water heat capacity" annotation (
        Dialog(tab="Dynamics", enable=not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

    // Variables
    Modelica.Units.SI.Efficiency eta=if effCur == Buildings.Fluid.Types.EfficiencyCurves.Constant
         then a[1] elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.Polynomial
         then Buildings.Utilities.Math.Functions.polynomial(a=a, x=y_internal)
         elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear
         then Buildings.Utilities.Math.Functions.quadraticLinear(
        a=aQuaLin,
        x1=y_internal,
        x2=MediumSte.saturationTemperature(port_a2.p)) else 0 "Boiler efficiency";
    Modelica.Units.SI.Power QFue_flow=y_internal*Q_flow_nominal/eta_nominal
      "Heat released by fuel";
    Modelica.Units.SI.Power QWat_flow=eta*QFue_flow
      "Heat transfer from gas into water";
    Modelica.Units.SI.MassFlowRate mFue_flow=QFue_flow/fue.h
      "Fuel mass flow rate";
    Modelica.Units.SI.VolumeFlowRate VFue_flow=mFue_flow/fue.d
      "Fuel volume flow rate";

    Modelica.Blocks.Interfaces.RealInput y(min=0, max=1)
      "Part load ratio"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealOutput VLiq(
      final quantity="Volume",
      final unit="m3",
      min=0) "Output liquid water volume"
      annotation (Placement(transformation(extent={{100,-90},{120,-70}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
                            "Heat port, can be used to connect to ambient"
      annotation (Placement(transformation(extent={{-10,90},{10,110}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heaCapDry(C=500*mDry,
        T(start=T_start))                       "Heat capacity of boiler metal"
      annotation (Placement(transformation(extent={{30,-78},{50,-58}})));
    Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation
      vol(
      redeclare final package MediumSte = MediumSte,
      redeclare final package MediumWat = MediumWat,
      final allowFlowReversal=allowFlowReversal2,
      final energyDynamics=energyDynamics,
      final massDynamics=massDynamics,
      final p_start=p_start,
      final m_flow_nominal=m2_flow_nominal,
      final show_T=show_T,
      final V=V,
      final fixed_p_start=fixed_p_start) "Steam/water control volume"
      annotation (Placement(transformation(extent={{8,-50},{28,-30}})));
    Buildings.Fluid.FixedResistances.PressureDrop res(
      redeclare final package Medium = MediumWat,
      final allowFlowReversal=allowFlowReversal2,
      final m_flow_nominal=m2_flow_nominal,
      final show_T=show_T,
      final dp_nominal=dp_nominal) "Flow resistance"
      annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));

    Modelica.Blocks.Interfaces.RealOutput QFueFlo(
      final quantity="HeatFlowRate",
      final unit="W",
      min=0) "Heat flow rate of the fuel"
      annotation (Placement(transformation(extent={{100,80},{120,100}})));

    Fluid.MixingVolumes.MixingVolume vol1(
      redeclare package Medium = Medium1,
      m_flow_nominal=m1_flow_nominal,
      V=V_com,
      nPorts=2) annotation (Placement(transformation(extent={{12,40},{-8,60}})));
    Fluid.Movers.FlowControlled_m_flow fan(
      redeclare package Medium = Medium1,
      m_flow_nominal=m1_flow_nominal,
      addPowerToMedium=false,
      nominalValuesDefineDefaultPressureCurve=true,
      tau=300,
      dp_nominal=1000)
      annotation (Placement(transformation(extent={{-60,50},{-40,30}})));

    Modelica.Blocks.Interfaces.RealOutput mFueFlo(
      final quantity="MassFlowRate",
      final unit="kg/s",
      min=0) "Mass flow rate of the fuel"
      annotation (Placement(transformation(extent={{100,66},{120,86}})));
  protected
    final parameter Boolean steadyDynamics=if energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState
         then true else false "= true, if steady state formulation";
    parameter Real eta_nominal(fixed=false)
      "Boiler efficiency at nominal condition";
    parameter Real aQuaLin[6]=if size(a, 1) == 6 then a else fill(0, 6)
      "Auxiliary variable for efficiency curve because quadraticLinear requires exactly 6 elements";

    Modelica.Blocks.Interfaces.RealInput y_internal(min=0, max=1)
      "Internal block needed for conditional input part load ratio";

    Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo
                            "Prescribed heat flow (if heatPort is connected)"
      annotation (Placement(transformation(extent={{-59,-90},{-39,-70}})));
    Modelica.Blocks.Sources.RealExpression Q_flow_in(y=QWat_flow)
      "Heat transfer from gas into water (if heatPort is connected)"
      annotation (Placement(transformation(extent={{-98,-90},{-78,-70}})));

    Modelica.Thermal.HeatTransfer.Components.ThermalConductor UAOve(G=UA)
      "Overall thermal conductance (if heatPort is connected)"
      annotation (Placement(transformation(extent={{-18,-74},{2,-54}})));

    Modelica.Blocks.Sources.RealExpression QFue_flow_out(y=QFue_flow)
      "Heat flow rate of the fuel"
      annotation (Placement(transformation(extent={{60,80},{80,100}})));
  public
    Modelica.Blocks.Sources.RealExpression Q_flow_exh(y=((1 - eta)*QFue_flow) +
          heatPort.Q_flow)  "Heat losses in the exhaust"
      annotation (Placement(transformation(extent={{80,50},{60,70}})));
  protected
    HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo1 if not steadyDynamics
      "Prescribed heat flow (if heatPort is connected)"
      annotation (Placement(transformation(extent={{41,50},{21,70}})));
  public
    Modelica.Blocks.Sources.RealExpression m_exh_flow(y=y*((mFueFlo*19.6) +
          mFueFlo))         "Heat losses in the exhaust"
      annotation (Placement(transformation(extent={{-100,-2},{-80,18}})));
  protected
    Modelica.Blocks.Sources.RealExpression mFue_flow_out(y=mFue_flow)
      "Mass flow rate of the fuel"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));
  initial equation
    if effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear then
      assert(size(a, 1) == 6, "The boiler has the efficiency curve set to 'Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear',
    and hence the parameter 'a' must have exactly 6 elements.
    However, only "   + String(size(a, 1)) + " elements were provided.");
    end if;

    if effCur == Buildings.Fluid.Types.EfficiencyCurves.Constant then
      eta_nominal = a[1];
    elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.Polynomial then
      eta_nominal = Buildings.Utilities.Math.Functions.polynomial(a=a, x=1);
    elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear then
      // For this efficiency curve, a must have 6 elements.
      eta_nominal = Buildings.Utilities.Math.Functions.quadraticLinear(
        a=aQuaLin,
        x1=1,
        x2=T_nominal);
    else
      eta_nominal = 999;
    end if;

  equation
    assert(eta > 0.001, "Efficiency curve is wrong.");

    connect(y, y_internal);

    if steadyDynamics then
      -QWat_flow = port_a2.m_flow*actualStream(port_a2.h_outflow) + port_b2.m_flow
        *actualStream(port_b2.h_outflow);
    end if;

    connect(UAOve.port_a, heatPort) annotation (Line(
        points={{-18,-64},{-20,-64},{-20,100},{0,100}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(Q_flow_in.y, preHeaFlo.Q_flow) annotation (Line(
        points={{-77,-80},{-59,-80}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(heaCapDry.port, UAOve.port_b) annotation (Line(points={{40,-78},{40,-80},
            {2,-80},{2,-64}}, color={191,0,0}));
    connect(preHeaFlo.port, UAOve.port_b)
      annotation (Line(points={{-39,-80},{2,-80},{2,-64}}, color={191,0,0}));
    connect(vol.heatPort, UAOve.port_b)
      annotation (Line(points={{18,-50},{18,-64},{2,-64}}, color={191,0,0}));
    connect(vol.port_b, port_b2)
      annotation (Line(points={{28,-40},{100,-40}}, color={0,127,255}));
    connect(port_a2, res.port_a)
      annotation (Line(points={{-100,-40},{-60,-40}}, color={0,127,255}));
    connect(res.port_b, vol.port_a)
      annotation (Line(points={{-40,-40},{8,-40}}, color={0,127,255}));
    connect(QFue_flow_out.y, QFueFlo)
      annotation (Line(points={{81,90},{110,90}}, color={0,0,127}));
    connect(vol1.ports[1], port_b1)
      annotation (Line(points={{3,40},{100,40}}, color={0,127,255}));
    connect(Q_flow_exh.y, preHeaFlo1.Q_flow) annotation (Line(
        points={{59,60},{41,60}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(preHeaFlo1.port, vol1.heatPort) annotation (Line(points={{21,60},{18,60},
            {18,50},{12,50}}, color={191,0,0}));
    connect(port_a1, fan.port_a)
      annotation (Line(points={{-100,40},{-60,40}}, color={0,127,255}));
    connect(fan.port_b, vol1.ports[2])
      annotation (Line(points={{-40,40},{1,40}}, color={0,127,255}));
    connect(vol.VLiq, VLiq) annotation (Line(points={{29,-33},{60,-33},{60,-80},{110,
            -80}}, color={0,0,127}));
    connect(m_exh_flow.y, fan.m_flow_in)
      annotation (Line(points={{-79,8},{-50,8},{-50,28}}, color={0,0,127}));
    connect(mFue_flow_out.y, mFueFlo)
      annotation (Line(points={{81,76},{110,76}}, color={0,0,127}));
    annotation (
      defaultComponentName="boi",
      Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-80,60},{80,-60}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-40,40},{40,-40}},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{20,18},{0,8},{20,-12},{0,-22}},
            color={0,0,0},
            smooth=Smooth.Bezier,
            extent={{-60,-22},{-36,2}}),
          Line(
            points={{-2,18},{-22,8},{-2,-12},{-22,-22}},
            color={0,0,0},
            smooth=Smooth.Bezier,
            extent={{-60,-22},{-36,2}})}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>
This model represents a steam boiler that discharges saturated 
steam and has an efficiency curve defined by a polynomial.
The efficiency in this model represents the fuel-to-water 
efficiency (e.g., thermal efficiency).
This model is similar to the 
<a href=\"modelica://Buildings.Fluid.Boilers.BoilerPolynomial\"> 
Buildings.Fluid.Boilers.BoilerPolynomial</a> for the efficiency 
and fuel mass flow rate computation with the following exceptions:
</p>
<ul>
<li>
Water enters <code>port_a</code> in liquid state and exits 
<code>port_b</code> in vapor state.
</li> 
<li>
The liquid and vapor phases are at equilibrium; thus, the steam
boiler is constrained to saturated states only with the volume
containing a wet steam mixture. 
</li>
<li>
If the boiler is configured in steady state, several blocks involving
the heat flow rate are conditionally removed to avoid overconstraining
the model. This is because the discharging fluid is constrained at 
a saturated state. The blocks that are conditionally removed as a 
result are within the green region in the below figure:
</li>
</ul>

<p align=\"center\">
<img src=\"modelica://Buildings/Resources/Images/Experimental/DHC/Plants/Steam/BaseClasses/BoilerPolynomial.png\" border=\"1\"
alt=\"Boiler polynomial steam with blocks in green conditionally removed if steady state\"/>
</p>
<h4>Implementation</h4>
<p>
In order to improve the numerical efficiency, this model follows 
the split-medium approach using the
<a href=\"modelica://Buildings.Fluid.Interfaces.PartialTwoPortTwoMedium\">
Buildings.Fluid.Interfaces.PartialTwoPortTwoMedium</a> interface model.
The saturated mixing volume for an evaporation process 
<a href=\"modelica://Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation\">
Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation</a> 
represents the phase change process of water from liquid 
to vapor at equilibrium.
</p>
<h4>Reference</h4>
<p>
Hinkelman, Kathryn, Saranya Anbarasu, Michael Wetter, 
Antoine Gautier, and Wangda Zuo. 2022. “A Fast and Accurate Modeling 
Approach for Water and Steam Thermodynamics with Practical 
Applications in District Heating System Simulation.” Preprint. February 24. 
<a href=\"http://dx.doi.org/10.13140/RG.2.2.20710.29762\">doi:10.13140/RG.2.2.20710.29762</a>.
</p>
</html>",   revisions="<html>
<ul>
<li>
February 25, 2022 by Kathryn Hinkelman:<br/>
Refactored base classes for improved extensibility and relocated models into Steam subpackages.
</li>
<li>
July 22, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"));
  end BoilerPolynomialExhaust;

  model BoilerPolynomialFurnaceHeatBalance "Variable efficiency with furnace heat balance function; A equilibrium boiler with water phase change from liquid to vapor, discharging
  saturated steam vapor, with the efficiency curve described by a polynomial."
    extends Buildings.BaseClasses.BaseIconLow;
    extends
      Buildings.Experimental.DHC.BaseClasses.Steam.PartialFourPortInterfaceThreeMedium(
      redeclare package Medium1 = MediumAir,
      redeclare package Medium2 = MediumWat,
      redeclare package Medium3 = MediumSte,
      final T_start=MediumSte.saturationTemperature(p_start));

    // Medium declarations
    replaceable package MediumAir =
        Buildings.Media.CombustionAir
      "Air medium on the combustion side";
    replaceable package MediumWat =
        Buildings.Media.Specialized.Water.TemperatureDependentDensity
      "Water medium - port_a2 (inlet)";
    replaceable package MediumSte = Buildings.Media.Steam
      "Steam medium - port_b2 (oulet)";
    // Initialization
    parameter Boolean fixed_p_start=false "Set to true if p_start is to be used as an explicit initial equation, 
    not an initial guess"   annotation (Dialog(tab="Initialization"));
    // Nominal conditions
    parameter Modelica.Units.SI.Power Q_flow_nominal "Nominal heating power";
    parameter Modelica.Units.SI.Temperature T_nominal=373.15 "Temperature used to compute nominal efficiency 
    (only used if efficiency curve depends on temperature)";

    // Efficiency, fuel, and boiler properties
    parameter Buildings.Fluid.Types.EfficiencyCurves effCur=Buildings.Fluid.Types.EfficiencyCurves.Constant
      "Curve used to compute the efficiency";
    parameter Real a[:]={0.9} "Coefficients for efficiency curve";
    parameter Buildings.Fluid.Data.Fuels.Generic fue "Fuel type"
      annotation (choicesAllMatching=true);
    parameter Modelica.Units.SI.ThermalConductance UA=0.05*Q_flow_nominal/30
      "Overall UA value";
    parameter Modelica.Units.SI.Volume V=1.5E-6*Q_flow_nominal
      "Total internal volume of boiler" annotation (Dialog(tab="Dynamics", enable=
           not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

    parameter Modelica.Units.SI.Volume V_com=m1_flow_nominal*tau1/rho1_nominal
      "Total internal volume of combustion" annotation (Dialog(tab="Dynamics",
          enable=not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

    parameter Modelica.Units.SI.Mass mDry=1.5E-3*Q_flow_nominal
      "Mass of boiler that will be lumped to water heat capacity" annotation (
        Dialog(tab="Dynamics", enable=not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

    // Variables
    Modelica.Units.SI.Efficiency eta=if effCur == Buildings.Fluid.Types.EfficiencyCurves.Constant
         then a[1] elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.Polynomial
         then Buildings.Utilities.Math.Functions.polynomial(a=a, x=y_internal)
         elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear
         then Buildings.Utilities.Math.Functions.quadraticLinear(
        a=aQuaLin,
        x1=y_internal,
        x2=MediumSte.saturationTemperature(port_a2.p)) else 0 "Boiler efficiency";
    Modelica.Units.SI.Power QFue_flow=y_internal*Q_flow_nominal/eta_nominal
      "Heat released by fuel";
    Modelica.Units.SI.Power QWat_flow=furHeaBal.eta_boi*QFue_flow
      "Heat transfer from gas into water";
    Modelica.Units.SI.MassFlowRate mFue_flow=QFue_flow/fue.h
      "Fuel mass flow rate";
    Modelica.Units.SI.VolumeFlowRate VFue_flow=mFue_flow/fue.d
      "Fuel volume flow rate";

    Modelica.Units.SI.MassFlowRate mExh_flow=mFue_flow+(mFue_flow/FA_ratio)
      "Fuel mass flow rate";

    Modelica.Blocks.Interfaces.RealInput y(min=0, max=1)
      "Part load ratio"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealOutput VLiq(
      final quantity="Volume",
      final unit="m3",
      min=0) "Output liquid water volume"
      annotation (Placement(transformation(extent={{100,-90},{120,-70}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
                            "Heat port, can be used to connect to ambient"
      annotation (Placement(transformation(extent={{-10,90},{10,110}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heaCapDry(C=500*mDry,
        T(start=T_start))                       "Heat capacity of boiler metal"
      annotation (Placement(transformation(extent={{40,-78},{60,-58}})));
    Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation
      vol(
      redeclare final package MediumSte = MediumSte,
      redeclare final package MediumWat = MediumWat,
      final allowFlowReversal=allowFlowReversal2,
      final energyDynamics=energyDynamics,
      final massDynamics=massDynamics,
      final p_start=p_start,
      final m_flow_nominal=m2_flow_nominal,
      final show_T=show_T,
      final V=V,
      final fixed_p_start=fixed_p_start) "Steam/water control volume"
      annotation (Placement(transformation(extent={{8,-50},{28,-30}})));
    Buildings.Fluid.FixedResistances.PressureDrop res(
      redeclare final package Medium = MediumWat,
      final allowFlowReversal=allowFlowReversal2,
      final m_flow_nominal=m2_flow_nominal,
      final show_T=show_T,
      final dp_nominal=dp_nominal) "Flow resistance"
      annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));

    Modelica.Blocks.Interfaces.RealOutput QFueFlo(
      final quantity="HeatFlowRate",
      final unit="W",
      min=0) "Heat flow rate of the fuel"
      annotation (Placement(transformation(extent={{100,80},{120,100}})));

    Fluid.MixingVolumes.MixingVolume vol1(
      redeclare package Medium = Medium1,
      energyDynamics=energyDynamics,
      massDynamics=massDynamics,
      mSenFac=mSenFac,
      m_flow_nominal=m1_flow_nominal,
      V=V_com,
      nPorts=2) annotation (Placement(transformation(extent={{12,40},{-8,60}})));
    Fluid.Movers.FlowControlled_m_flow fan(
      redeclare package Medium = Medium1,
      m_flow_nominal=m_flow_nominal,
      addPowerToMedium=false,
      nominalValuesDefineDefaultPressureCurve=true,
      tau=30,
      dp_nominal=dp_nominal)
      annotation (Placement(transformation(extent={{-80,50},{-60,30}})));

    Modelica.Blocks.Interfaces.RealOutput mFueFlo(
      final quantity="MassFlowRate",
      final unit="kg/s",
      min=0) "Mass flow rate of the fuel"
      annotation (Placement(transformation(extent={{100,66},{120,86}})));
  protected
    final parameter Boolean steadyDynamics=if energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState
         then true else false "= true, if steady state formulation";
    parameter Real eta_nominal(fixed=false)
      "Boiler efficiency at nominal condition";
    parameter Real aQuaLin[6]=if size(a, 1) == 6 then a else fill(0, 6)
      "Auxiliary variable for efficiency curve because quadraticLinear requires exactly 6 elements";

    Modelica.Blocks.Interfaces.RealInput y_internal(min=0, max=1)
      "Internal block needed for conditional input part load ratio";

    Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo
                            "Prescribed heat flow (if heatPort is connected)"
      annotation (Placement(transformation(extent={{-59,-90},{-39,-70}})));
  public
    Modelica.Blocks.Sources.RealExpression Q_Wat_flow(y=QWat_flow)
      "Heat transfer from gas into water (if heatPort is connected)"
      annotation (Placement(transformation(extent={{-100,-90},{-80,-70}})));

  protected
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor UAOve(G=UA)
      "Overall thermal conductance (if heatPort is connected)"
      annotation (Placement(transformation(extent={{-20,-70},{0,-50}})));

    Modelica.Blocks.Sources.RealExpression QFue_flow_out(y=QFue_flow)
      "Heat flow rate of the fuel"
      annotation (Placement(transformation(extent={{60,80},{80,100}})));
  public
    Modelica.Blocks.Sources.RealExpression qExhFlo(y=(furHeaBal.per_exh*
          QFue_flow)) "Heat losses in the exhaust"
      annotation (Placement(transformation(extent={{80,40},{60,60}})));
  protected
    HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo1 if not steadyDynamics
      "Prescribed heat flow (if heatPort is connected)"
      annotation (Placement(transformation(extent={{41,40},{21,60}})));
  public
    Modelica.Blocks.Sources.RealExpression mExhFlow(y=0.9*(((mFueFlo/
          furHeaBal.fs) + mFueFlo)))
                       "Heat losses in the exhaust"
      annotation (Placement(transformation(extent={{-100,-2},{-80,18}})));
    Fluid.Sensors.TemperatureTwoPort
                                   senTem(
      redeclare package Medium = Medium1,
        m_flow_nominal=m1_flow_nominal,
      T_start(displayUnit="K"))
      annotation (Placement(transformation(extent={{-50,50},{-30,30}})));
    BaseClasses.FurnaceHeatBalance furHeaBal(
      T_exh_nominal=T_exh_nominal,
      QFue=fue.h,
      FA_ratio=FA_ratio)
      annotation (Placement(transformation(extent={{0,12},{20,-8}})));
  public
    Modelica.Blocks.Sources.RealExpression qLos(y=heatPort.Q_flow)
      "Heat losses ifrom boiler casing"
      annotation (Placement(transformation(extent={{-60,-14},{-40,6}})));
    parameter Real FA_ratio
      "Fuel air ratio, alpha (20% excess air, FA_ratio = 1.20)";
    parameter Modelica.Units.SI.Temperature T_exh_nominal=373.15
      "Exhaust temperature used to compute nominal efficiency";
    parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=8
      "Nominal mass flow rate" annotation (Dialog(group="Compressor Fan"));
    parameter Modelica.Units.SI.PressureDifference dp_nominal=1000
      "Nominal pressure raise, used for default pressure curve if not specified in record per"
      annotation (Dialog(group="Compressor Fan"));
    parameter Real mSenFac=1
      "Factor for scaling the sensible thermal mass of the volume";
  protected
    Modelica.Blocks.Sources.RealExpression mFue_flow_out(y=mFue_flow)
      "Mass flow rate of the fuel"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));
  initial equation
    if effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear then
      assert(size(a, 1) == 6, "The boiler has the efficiency curve set to 'Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear',
    and hence the parameter 'a' must have exactly 6 elements.
    However, only "   + String(size(a, 1)) + " elements were provided.");
    end if;

    if effCur == Buildings.Fluid.Types.EfficiencyCurves.Constant then
      eta_nominal = a[1];
    elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.Polynomial then
      eta_nominal = Buildings.Utilities.Math.Functions.polynomial(a=a, x=1);
    elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear then
      // For this efficiency curve, a must have 6 elements.
      eta_nominal = Buildings.Utilities.Math.Functions.quadraticLinear(
        a=aQuaLin,
        x1=1,
        x2=T_nominal);
    else
      eta_nominal = 999;
    end if;

  equation
    assert(eta > 0.001, "Efficiency curve is wrong.");

    connect(y, y_internal);

    if steadyDynamics then
      -QWat_flow = port_a2.m_flow*actualStream(port_a2.h_outflow) + port_b2.m_flow
        *actualStream(port_b2.h_outflow);
    end if;

    connect(UAOve.port_a, heatPort) annotation (Line(
        points={{-20,-60},{-20,80},{0,80},{0,100}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(Q_Wat_flow.y, preHeaFlo.Q_flow) annotation (Line(
        points={{-79,-80},{-59,-80}},
        color={0,0,127},
        smooth=Smooth.None,
        pattern=LinePattern.Dash));
    connect(heaCapDry.port, UAOve.port_b) annotation (Line(points={{50,-78},{50,
            -80},{18,-80},{18,-60},{0,-60}},
                              color={191,0,0}));
    connect(preHeaFlo.port, UAOve.port_b)
      annotation (Line(points={{-39,-80},{18,-80},{18,-60},{0,-60}},
                                                           color={191,0,0}));
    connect(vol.heatPort, UAOve.port_b)
      annotation (Line(points={{18,-50},{18,-60},{0,-60}}, color={191,0,0}));
    connect(vol.port_b, port_b2)
      annotation (Line(points={{28,-40},{100,-40}}, color={0,127,255},
        thickness=0.5));
    connect(port_a2, res.port_a)
      annotation (Line(points={{-100,-40},{-60,-40}}, color={0,127,255},
        thickness=0.5));
    connect(res.port_b, vol.port_a)
      annotation (Line(points={{-40,-40},{8,-40}}, color={0,127,255},
        thickness=0.5));
    connect(QFue_flow_out.y, QFueFlo)
      annotation (Line(points={{81,90},{110,90}}, color={0,0,127},
        pattern=LinePattern.Dash));
    connect(vol1.ports[1], port_b1)
      annotation (Line(points={{3,40},{100,40}}, color={0,127,255},
        thickness=0.5));
    connect(qExhFlo.y, preHeaFlo1.Q_flow) annotation (Line(
        points={{59,50},{41,50}},
        color={0,0,127},
        smooth=Smooth.None,
        pattern=LinePattern.Dash));
    connect(preHeaFlo1.port, vol1.heatPort) annotation (Line(points={{21,50},{12,50}},
                              color={191,0,0}));
    connect(port_a1, fan.port_a)
      annotation (Line(points={{-100,40},{-80,40}}, color={0,127,255},
        thickness=0.5));
    connect(vol.VLiq, VLiq) annotation (Line(points={{29,-33},{80,-33},{80,-80},
            {110,-80}},
                   color={0,0,127},
        pattern=LinePattern.Dash));
    connect(mExhFlow.y, fan.m_flow_in) annotation (Line(
        points={{-79,8},{-70,8},{-70,28}},
        color={0,0,127},
        pattern=LinePattern.Dash));
    connect(mFue_flow_out.y, mFueFlo)
      annotation (Line(points={{81,76},{110,76}}, color={0,0,127},
        pattern=LinePattern.Dash));
    connect(fan.port_b, senTem.port_a)
      annotation (Line(points={{-60,40},{-50,40}}, color={0,127,255},
        thickness=0.5));
    connect(senTem.port_b, vol1.ports[2])
      annotation (Line(points={{-30,40},{1,40}}, color={0,127,255},
        thickness=0.5));
    connect(senTem.T, furHeaBal.T_Air_in) annotation (Line(
        points={{-40,29},{-40,4},{-2,4}},
        color={0,0,127},
        pattern=LinePattern.Dash));
    connect(qLos.y, furHeaBal.qLos) annotation (Line(
        points={{-39,-4},{-2,-4}},
        color={0,0,127},
        pattern=LinePattern.Dash));
    annotation (
      defaultComponentName="boi",
      Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-80,60},{80,-60}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-40,40},{40,-40}},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{20,18},{0,8},{20,-12},{0,-22}},
            color={0,0,0},
            smooth=Smooth.Bezier,
            extent={{-60,-22},{-36,2}}),
          Line(
            points={{-2,18},{-22,8},{-2,-12},{-22,-22}},
            color={0,0,0},
            smooth=Smooth.Bezier,
            extent={{-60,-22},{-36,2}})}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>
This model represents a steam boiler that discharges saturated 
steam and has an efficiency curve defined by a polynomial.
The efficiency in this model represents the fuel-to-water 
efficiency (e.g., thermal efficiency).
This model is similar to the 
<a href=\"modelica://Buildings.Fluid.Boilers.BoilerPolynomial\"> 
Buildings.Fluid.Boilers.BoilerPolynomial</a> for the efficiency 
and fuel mass flow rate computation with the following exceptions:
</p>
<ul>
<li>
Water enters <code>port_a</code> in liquid state and exits 
<code>port_b</code> in vapor state.
</li> 
<li>
The liquid and vapor phases are at equilibrium; thus, the steam
boiler is constrained to saturated states only with the volume
containing a wet steam mixture. 
</li>
<li>
If the boiler is configured in steady state, several blocks involving
the heat flow rate are conditionally removed to avoid overconstraining
the model. This is because the discharging fluid is constrained at 
a saturated state. The blocks that are conditionally removed as a 
result are within the green region in the below figure:
</li>
</ul>

<p align=\"center\">
<img src=\"modelica://Buildings/Resources/Images/Experimental/DHC/Plants/Steam/BaseClasses/BoilerPolynomial.png\" border=\"1\"
alt=\"Boiler polynomial steam with blocks in green conditionally removed if steady state\"/>
</p>
<h4>Implementation</h4>
<p>
In order to improve the numerical efficiency, this model follows 
the split-medium approach using the
<a href=\"modelica://Buildings.Fluid.Interfaces.PartialTwoPortTwoMedium\">
Buildings.Fluid.Interfaces.PartialTwoPortTwoMedium</a> interface model.
The saturated mixing volume for an evaporation process 
<a href=\"modelica://Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation\">
Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation</a> 
represents the phase change process of water from liquid 
to vapor at equilibrium.
</p>
<h4>Reference</h4>
<p>
Hinkelman, Kathryn, Saranya Anbarasu, Michael Wetter, 
Antoine Gautier, and Wangda Zuo. 2022. “A Fast and Accurate Modeling 
Approach for Water and Steam Thermodynamics with Practical 
Applications in District Heating System Simulation.” Preprint. February 24. 
<a href=\"http://dx.doi.org/10.13140/RG.2.2.20710.29762\">doi:10.13140/RG.2.2.20710.29762</a>.
</p>
</html>",   revisions="<html>
<ul>
<li>
February 25, 2022 by Kathryn Hinkelman:<br/>
Refactored base classes for improved extensibility and relocated models into Steam subpackages.
</li>
<li>
July 22, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"));
  end BoilerPolynomialFurnaceHeatBalance;

  package Validation
    extends Modelica.Icons.ExamplesPackage;

    package Boiler
      model BoilerPolynomialFurnaceHeatBalance
        "Example model for the steam boiler with a polynomial efficiency curve"
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
        parameter Modelica.Units.SI.Power Q_flow_nominal = 29306000 "Nominal power";
        parameter Modelica.Units.SI.SpecificEnthalpy dh_nominal=
          MediumSte.specificEnthalpy(
            MediumSte.setState_pTX(p=p_nominal, T=T_nominal, X=MediumSte.X_default))
          "Nominal change in enthalpy";
        parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=
          Q_flow_nominal/dh_nominal/2
          "Nominal mass flow rate";
        parameter Modelica.Units.SI.PressureDifference dp_nominal = 185787
          "Pressure drop at m_flow_nominal";

        Buildings.Fluid.Sources.Boundary_pT sin(
          redeclare package Medium = MediumSte,
          p(displayUnit="bar") = 903213,
          T(displayUnit="K") = 427,
          nPorts=1)
          "Sink"
          annotation (Placement(transformation(extent={{80,10},{60,30}})));
        Buildings.Fluid.Sources.Boundary_pT sou(
          redeclare package Medium = MediumWat,
          p(displayUnit="bar") = 1758000,
          T(displayUnit="K") = 375.9,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
        Buildings.HeatTransfer.Sources.FixedTemperature TAmb(T=303.15)
          "Ambient temperature in boiler room"
          annotation (Placement(transformation(extent={{-30,102},{-10,122}})));
        .Buildings.GEDHeatingSystemCalibration.CUBoulder.Components.BoilerPolynomialFurnaceHeatBalance
          boiDyn(
          m1_flow_nominal=7.5,
          m2_flow_nominal=6.8,
          redeclare package MediumAir = MediumAir,
          redeclare package MediumSte = MediumSte,
          redeclare package MediumWat = MediumWat,
          energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
          Q_flow_nominal=Q_flow_nominal,
          a={0.8},
          fue=Buildings.Fluid.Data.Fuels.Generic(
                        h=47522000,
                        d=0.8,
                        mCO2=2),
          dp_nominal=dp_nominal,
          UA=0.05*Q_flow_nominal/100,
          V=25,
          V_com=21.94,
          FA_ratio=1.15,
          T_exh_nominal(displayUnit="K") = 423)
          "Steam boiler with dynamic balance"
          annotation (Placement(transformation(extent={{-10,42},{10,62}})));
        Fluid.Sources.Boundary_pT           sou1(
          redeclare package Medium = MediumAir,
          use_T_in=true,
          T=573.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-80,70},{-60,90}})));
        Fluid.Sources.Boundary_pT           sou2(redeclare package Medium = MediumAir,
          T(displayUnit="K") = 577,
            nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{80,70},{60,90}})));

        Modelica.Blocks.Sources.Constant
                                     const(k=309)
          annotation (Placement(transformation(extent={{-120,74},{-100,94}})));
        Fluid.Movers.FlowControlled_m_flow fwPum(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m_flow_nominal,
          addPowerToMedium=false,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal=dp_nominal) "Feed water pump"
          annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
        Modelica.Blocks.Sources.Constant
                                     const1(k=5.95)
          annotation (Placement(transformation(extent={{-120,30},{-100,50}})));
        Modelica.Blocks.Sources.SawTooth
                                      sawTooth(
          amplitude=1,
          period=500)
          annotation (Placement(transformation(extent={{-120,114},{-100,134}})));
        Modelica.Blocks.Sources.RealExpression QFue(y=boiDyn.QFue_flow)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{-170,-66},{-150,-46}})));
        Modelica.Blocks.Sources.RealExpression QWat(y=boiDyn.QWat_flow)
          "Water heat flow rate"
          annotation (Placement(transformation(extent={{-170,-116},{-150,-96}})));
        Modelica.Blocks.Sources.RealExpression QLoss(y=boiDyn.heatPort.Q_flow)
          "Loss from boiler casing"
          annotation (Placement(transformation(extent={{-62,-118},{-42,-98}})));
        Modelica.Blocks.Sources.RealExpression QFlue(y=boiDyn.qExhFlo.y)
          "Heat losses into flue gas"
          annotation (Placement(transformation(extent={{-62,-68},{-42,-48}})));
        Modelica.Blocks.Sources.RealExpression Qcomb(y=((boiDyn.port_b1.h_outflow -
              boiDyn.port_a1.h_outflow)*boiDyn.port_a1.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{66,-68},{86,-48}})));
        Modelica.Blocks.Sources.RealExpression Qevap(y=((boiDyn.port_b2.h_outflow -
              boiDyn.port_a2.h_outflow)*boiDyn.port_b2.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{66,-118},{86,-98}})));
        Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-132,-66},{-112,-46}})));
        Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-132,-116},{-112,-96}})));
        Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-20,-68},{0,-48}})));
        Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-20,-118},{0,-98}})));
        Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{114,-68},{134,-48}})));
        Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{114,-118},{134,-98}})));
        Modelica.Blocks.Math.Division perQWat
          annotation (Placement(transformation(extent={{-96,-90},{-76,-110}})));
        Modelica.Blocks.Math.Division perQLoss
          annotation (Placement(transformation(extent={{14,-92},{34,-112}})));
        Modelica.Blocks.Math.Division perQFlue
          annotation (Placement(transformation(extent={{16,-74},{36,-54}})));
        Modelica.Blocks.Math.Division perQcomb
          annotation (Placement(transformation(extent={{150,-74},{170,-54}})));
        Modelica.Blocks.Math.Division perQeva
          annotation (Placement(transformation(extent={{150,-92},{170,-112}})));
        Modelica.Blocks.Math.Add diffQevap(k1=-1)
          annotation (Placement(transformation(extent={{32,-146},{12,-166}})));
        Modelica.Blocks.Math.Add diffQcomb(k1=-1)
          annotation (Placement(transformation(extent={{90,-146},{70,-166}})));
      equation
        connect(TAmb.port, boiDyn.heatPort)
          annotation (Line(points={{-10,112},{0,112},{0,62}},   color={191,0,0}));
        connect(boiDyn.port_b2, sin.ports[1]) annotation (Line(points={{10,48},{20,48},
                {20,20},{60,20}},   color={0,127,255}));
        connect(sou1.ports[1], boiDyn.port_a1) annotation (Line(points={{-60,80},{-42,
                80},{-42,56},{-10,56}},
                                    color={0,127,255}));
        connect(boiDyn.port_b1, sou2.ports[1]) annotation (Line(points={{10,56},{20,
                56},{20,80},{60,80}},
                                  color={0,127,255}));
        connect(const.y, sou1.T_in)
          annotation (Line(points={{-99,84},{-82,84}}, color={0,0,127}));
        connect(sou.ports[1], fwPum.port_a)
          annotation (Line(points={{-60,20},{-40,20}},   color={0,127,255}));
        connect(fwPum.port_b, boiDyn.port_a2) annotation (Line(points={{-20,20},{-16,
                20},{-16,48},{-10,48}},  color={0,127,255}));
        connect(const1.y, fwPum.m_flow_in)
          annotation (Line(points={{-99,40},{-30,40},{-30,32}},    color={0,0,127}));
        connect(sawTooth.y, boiDyn.y) annotation (Line(points={{-99,124},{-40,
                124},{-40,60},{-12,60}}, color={0,0,127}));
        connect(QFue.y, iQfue.u)
          annotation (Line(points={{-149,-56},{-134,-56}}, color={0,0,127}));
        connect(QWat.y, iQWat.u)
          annotation (Line(points={{-149,-106},{-134,-106}}, color={0,0,127}));
        connect(QFlue.y, iQFlue.u)
          annotation (Line(points={{-41,-58},{-22,-58}}, color={0,0,127}));
        connect(QLoss.y, iQLoss.u)
          annotation (Line(points={{-41,-108},{-22,-108}}, color={0,0,127}));
        connect(iQComb.u, Qcomb.y)
          annotation (Line(points={{112,-58},{87,-58}}, color={0,0,127}));
        connect(iQevap.u, Qevap.y)
          annotation (Line(points={{112,-108},{87,-108}}, color={0,0,127}));
        connect(iQWat.y, perQWat.u1)
          annotation (Line(points={{-111,-106},{-98,-106}}, color={0,0,127}));
        connect(iQfue.y, perQWat.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-94},{-98,-94}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQLoss.y, perQLoss.u1)
          annotation (Line(points={{1,-108},{12,-108}}, color={0,0,127}));
        connect(iQfue.y, perQLoss.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{6,-82},{6,-96},{12,-96}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQFlue.y, perQFlue.u1)
          annotation (Line(points={{1,-58},{14,-58}}, color={0,0,127}));
        connect(iQfue.y, perQFlue.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{6,-82},{6,-70},{14,-70}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQComb.y, perQcomb.u1)
          annotation (Line(points={{135,-58},{148,-58}}, color={0,0,127}));
        connect(iQfue.y, perQcomb.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{140,-82},{140,-70},{148,-70}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQevap.y, perQeva.u1)
          annotation (Line(points={{135,-108},{148,-108}}, color={0,0,127}));
        connect(iQfue.y, perQeva.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{140,-82},{140,-96},{148,-96}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(perQcomb.y, diffQcomb.u2) annotation (Line(points={{171,-64},{178,-64},
                {178,-150},{92,-150}}, color={0,0,127}));
        connect(perQeva.y, diffQevap.u2) annotation (Line(points={{171,-102},{170,-102},
                {170,-142},{34,-142},{34,-150}}, color={0,0,127}));
        connect(perQWat.y, diffQevap.u1) annotation (Line(points={{-75,-100},{-68,-100},
                {-68,-176},{46,-176},{46,-162},{34,-162}}, color={244,125,35}));
        connect(perQFlue.y, diffQcomb.u1) annotation (Line(points={{37,-64},{60,-64},{
                60,-176},{102,-176},{102,-162},{92,-162}}, color={244,125,35}));
        annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Plants/Steam/BaseClasses/Examples/BoilerPolynomialFourPort.mos"
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
</html>",       revisions="<html>
<ul>
<li>
July 23, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"),Diagram(coordinateSystem(extent={{-180,-100},{180,160}}), graphics={
                Rectangle(
                extent={{-180,-40},{180,-180}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid), Text(
                extent={{-176,-26},{-70,-46}},
                textColor={28,108,200},
                textString="Steady state verification (Heat balance)")}));
      end BoilerPolynomialFurnaceHeatBalance;

      model FurnaceHeatBalance
          extends Modelica.Icons.Example;
        Buildings.GEDHeatingSystemCalibration.CUBoulder.Components.BaseClasses.FurnaceHeatBalance
          furnaceHeatBalance(
          T_exh_nominal(displayUnit="K") = 577,
          QFue=47522000,
          FA_ratio=1.15)
          annotation (Placement(transformation(extent={{-20,0},{0,20}})));
        Modelica.Blocks.Sources.Ramp     ramp(
          height=22,
          duration=3600,
          offset=277.5)
          annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
        Modelica.Blocks.Sources.Constant const(k=2742000)
          annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      equation

        connect(const.y, furnaceHeatBalance.qLos) annotation (Line(points={{-59,30},{-28,
                30},{-28,16},{-22,16}}, color={0,0,127}));
        connect(ramp.y, furnaceHeatBalance.T_Air_in) annotation (Line(points={{
                -59,-10},{-28,-10},{-28,8},{-22,8}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end FurnaceHeatBalance;

      model BoilerPolynomialExhaust
        "Example model for the steam boiler with a polynomial efficiency curve"
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
        parameter Modelica.Units.SI.Power Q_flow_nominal = 29306000 "Nominal power";
        parameter Modelica.Units.SI.SpecificEnthalpy dh_nominal=
          MediumSte.specificEnthalpy(
            MediumSte.setState_pTX(p=p_nominal, T=T_nominal, X=MediumSte.X_default))
          "Nominal change in enthalpy";
        parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=
          Q_flow_nominal/dh_nominal/2
          "Nominal mass flow rate";
        parameter Modelica.Units.SI.PressureDifference dp_nominal = 800000
          "Pressure drop at m_flow_nominal";

        Buildings.Fluid.Sources.Boundary_pT sin(
          redeclare package Medium = MediumSte,
          p(displayUnit="bar") = 300000,
          T=423.15,
          nPorts=1)
          "Sink"
          annotation (Placement(transformation(extent={{80,10},{60,30}})));
        Buildings.Fluid.Sources.Boundary_pT sou(
          redeclare package Medium = MediumWat,
          use_T_in=true,
          T=303.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
        Buildings.HeatTransfer.Sources.FixedTemperature TAmb(T=288.15)
          "Ambient temperature in boiler room"
          annotation (Placement(transformation(extent={{-30,102},{-10,122}})));
        GEDHeatingSystemCalibration.CUBoulder.Components.BoilerPolynomialExhaust
          boiDyn(
          m1_flow_nominal=1,
          m2_flow_nominal=m_flow_nominal,
          redeclare package MediumAir = MediumAir,
          redeclare package MediumSte = MediumSte,
          redeclare package MediumWat = MediumWat,
          energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
          Q_flow_nominal=Q_flow_nominal,
          a={0.8},
          fue=Buildings.Fluid.Data.Fuels.NaturalGasHigherHeatingValue(),
          dp_nominal=dp_nominal,
          UA=0.05*Q_flow_nominal/100,
          V=25,
          V_com=21.94) "Steam boiler with dynamic balance"
          annotation (Placement(transformation(extent={{-10,42},{10,62}})));
        Fluid.Sources.Boundary_pT           sou1(
          redeclare package Medium = MediumAir,
          use_T_in=true,
          T=573.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-80,70},{-60,90}})));
        Fluid.Sources.Boundary_pT           sou2(redeclare package Medium = MediumAir,
            nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{80,70},{60,90}})));
        Modelica.Blocks.Sources.Ramp ramp(
          height=300,
          duration=3000,
          offset=273.15 + 25)
          annotation (Placement(transformation(extent={{-120,74},{-100,94}})));
        Fluid.Movers.FlowControlled_m_flow fwPum(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m_flow_nominal,
          addPowerToMedium=false,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal=dp_nominal) "Feed water pump"
          annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
        Modelica.Blocks.Sources.Constant
                                     const1(k=13.3)
          annotation (Placement(transformation(extent={{-120,30},{-100,50}})));
        Modelica.Blocks.Sources.Constant
                                      const(k=1)
          annotation (Placement(transformation(extent={{-120,114},{-100,134}})));
        Modelica.Blocks.Sources.RealExpression QFue(y=boiDyn.QFue_flow)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{-170,-66},{-150,-46}})));
        Modelica.Blocks.Sources.RealExpression QWat(y=boiDyn.QWat_flow)
          "Water heat flow rate"
          annotation (Placement(transformation(extent={{-170,-116},{-150,-96}})));
        Modelica.Blocks.Sources.RealExpression QLoss(y=boiDyn.heatPort.Q_flow)
          "Loss from boiler casing"
          annotation (Placement(transformation(extent={{-62,-118},{-42,-98}})));
        Modelica.Blocks.Sources.RealExpression QFlue(y=boiDyn.Q_flow_exh.y)
          "Heat losses into flue gas"
          annotation (Placement(transformation(extent={{-62,-68},{-42,-48}})));
        Modelica.Blocks.Sources.RealExpression Qcomb(y=((boiDyn.port_b1.h_outflow -
              boiDyn.port_a1.h_outflow)*boiDyn.port_a1.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{66,-68},{86,-48}})));
        Modelica.Blocks.Sources.RealExpression Qevap(y=((boiDyn.port_b2.h_outflow -
              boiDyn.port_a2.h_outflow)*boiDyn.port_a2.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{66,-118},{86,-98}})));
        Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-132,-66},{-112,-46}})));
        Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-132,-116},{-112,-96}})));
        Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-20,-68},{0,-48}})));
        Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-20,-118},{0,-98}})));
        Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{114,-68},{134,-48}})));
        Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{114,-118},{134,-98}})));
        Modelica.Blocks.Math.Division perQWat
          annotation (Placement(transformation(extent={{-96,-90},{-76,-110}})));
        Modelica.Blocks.Math.Division perQLoss
          annotation (Placement(transformation(extent={{14,-92},{34,-112}})));
        Modelica.Blocks.Math.Division perQFlue
          annotation (Placement(transformation(extent={{16,-74},{36,-54}})));
        Modelica.Blocks.Math.Division perQcomb
          annotation (Placement(transformation(extent={{150,-74},{170,-54}})));
        Modelica.Blocks.Math.Division perQeva
          annotation (Placement(transformation(extent={{150,-92},{170,-112}})));
        Modelica.Blocks.Math.Add diffQevap(k1=-1)
          annotation (Placement(transformation(extent={{32,-146},{12,-166}})));
        Modelica.Blocks.Math.Add diffQcomb(k1=-1)
          annotation (Placement(transformation(extent={{90,-146},{70,-166}})));
        Modelica.Blocks.Sources.Ramp ramp1(
          height=100,
          duration=3000,
          offset=273.15 + 25)
          annotation (Placement(transformation(extent={{-156,-8},{-136,12}})));
      equation
        connect(TAmb.port, boiDyn.heatPort)
          annotation (Line(points={{-10,112},{0,112},{0,62}},   color={191,0,0}));
        connect(boiDyn.port_b2, sin.ports[1]) annotation (Line(points={{10,48},{20,48},
                {20,20},{60,20}},   color={0,127,255}));
        connect(sou1.ports[1], boiDyn.port_a1) annotation (Line(points={{-60,80},{-42,
                80},{-42,56},{-10,56}},
                                    color={0,127,255}));
        connect(boiDyn.port_b1, sou2.ports[1]) annotation (Line(points={{10,56},{20,
                56},{20,80},{60,80}},
                                  color={0,127,255}));
        connect(ramp.y, sou1.T_in)
          annotation (Line(points={{-99,84},{-82,84}}, color={0,0,127}));
        connect(sou.ports[1], fwPum.port_a)
          annotation (Line(points={{-60,20},{-40,20}},   color={0,127,255}));
        connect(fwPum.port_b, boiDyn.port_a2) annotation (Line(points={{-20,20},{-16,
                20},{-16,48},{-10,48}},  color={0,127,255}));
        connect(const1.y, fwPum.m_flow_in)
          annotation (Line(points={{-99,40},{-30,40},{-30,32}},    color={0,0,127}));
        connect(const.y, boiDyn.y) annotation (Line(points={{-99,124},{-40,124},{-40,
                60},{-12,60}}, color={0,0,127}));
        connect(QFue.y, iQfue.u)
          annotation (Line(points={{-149,-56},{-134,-56}}, color={0,0,127}));
        connect(QWat.y, iQWat.u)
          annotation (Line(points={{-149,-106},{-134,-106}}, color={0,0,127}));
        connect(QFlue.y, iQFlue.u)
          annotation (Line(points={{-41,-58},{-22,-58}}, color={0,0,127}));
        connect(QLoss.y, iQLoss.u)
          annotation (Line(points={{-41,-108},{-22,-108}}, color={0,0,127}));
        connect(iQComb.u, Qcomb.y)
          annotation (Line(points={{112,-58},{87,-58}}, color={0,0,127}));
        connect(iQevap.u, Qevap.y)
          annotation (Line(points={{112,-108},{87,-108}}, color={0,0,127}));
        connect(iQWat.y, perQWat.u1)
          annotation (Line(points={{-111,-106},{-98,-106}}, color={0,0,127}));
        connect(iQfue.y, perQWat.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-94},{-98,-94}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQLoss.y, perQLoss.u1)
          annotation (Line(points={{1,-108},{12,-108}}, color={0,0,127}));
        connect(iQfue.y, perQLoss.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{6,-82},{6,-96},{12,-96}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQFlue.y, perQFlue.u1)
          annotation (Line(points={{1,-58},{14,-58}}, color={0,0,127}));
        connect(iQfue.y, perQFlue.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{6,-82},{6,-70},{14,-70}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQComb.y, perQcomb.u1)
          annotation (Line(points={{135,-58},{148,-58}}, color={0,0,127}));
        connect(iQfue.y, perQcomb.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{140,-82},{140,-70},{148,-70}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQevap.y, perQeva.u1)
          annotation (Line(points={{135,-108},{148,-108}}, color={0,0,127}));
        connect(iQfue.y, perQeva.u2) annotation (Line(
            points={{-111,-56},{-104,-56},{-104,-82},{140,-82},{140,-96},{148,-96}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(perQcomb.y, diffQcomb.u2) annotation (Line(points={{171,-64},{178,-64},
                {178,-150},{92,-150}}, color={0,0,127}));
        connect(perQeva.y, diffQevap.u2) annotation (Line(points={{171,-102},{170,-102},
                {170,-142},{34,-142},{34,-150}}, color={0,0,127}));
        connect(perQWat.y, diffQevap.u1) annotation (Line(points={{-75,-100},{-68,-100},
                {-68,-176},{46,-176},{46,-162},{34,-162}}, color={244,125,35}));
        connect(perQFlue.y, diffQcomb.u1) annotation (Line(points={{37,-64},{60,-64},{
                60,-176},{102,-176},{102,-162},{92,-162}}, color={244,125,35}));
        connect(ramp1.y, sou.T_in) annotation (Line(points={{-135,2},{-90,2},{-90,24},
                {-82,24}}, color={0,0,127}));
        annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Plants/Steam/BaseClasses/Examples/BoilerPolynomialFourPort.mos"
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
</html>",       revisions="<html>
<ul>
<li>
July 23, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"),Diagram(coordinateSystem(extent={{-180,-100},{180,140}}), graphics={
                Rectangle(
                extent={{-180,-40},{180,-180}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid), Text(
                extent={{-176,-26},{-70,-46}},
                textColor={28,108,200},
                textString="Steady state verification (Heat balance)")}));
      end BoilerPolynomialExhaust;
    end Boiler;

    model HeatBalanceBoiler
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

      Fluid.Sources.Boundary_pT           sou(
        redeclare package Medium = MediumWat,
        use_T_in=true,
        T=303.15,
        nPorts=1)
        "Source"
        annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
      Fluid.Movers.FlowControlled_m_flow fwPum(
        redeclare package Medium = MediumWat,
        m_flow_nominal=m_flow_nominal,
        addPowerToMedium=false,
        nominalValuesDefineDefaultPressureCurve=true,
        dp_nominal=dp_nominal) "Feed water pump"
        annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
      Modelica.Blocks.Sources.Constant
                                   mFloFW(k=4)
        annotation (Placement(transformation(extent={{-120,-10},{-100,10}})));
      Modelica.Blocks.Sources.Constant
                                   Tfw_in(k=300)
        annotation (Placement(transformation(extent={{-118,-50},{-98,-30}})));
      Modelica.Blocks.Sources.Ramp     Tair_in(
        height=100,
        duration=3600,
        offset=290)
        annotation (Placement(transformation(extent={{-120,34},{-100,54}})));
      Fluid.Sources.Boundary_pT pro(
        redeclare package Medium = MediumAir,
        use_T_in=true,
        T=573.15,
        nPorts=1) "Source"
        annotation (Placement(transformation(extent={{-80,30},{-60,50}})));
      HeatTransfer.Sources.FixedTemperature           TAmb(T(displayUnit="K") = 284)
        "Ambient temperature in boiler room"
        annotation (Placement(transformation(extent={{-14,84},{6,104}})));
      Fluid.Sources.Boundary_pT exh(
        redeclare package Medium = MediumAir,
        T(displayUnit="K"),
        nPorts=1) "Source"
        annotation (Placement(transformation(extent={{100,30},{80,50}})));

      Fluid.Sources.Boundary_pT           sin(
        redeclare package Medium = MediumSte,
        p(displayUnit="bar") = 900000,
        T=453.15,
        nPorts=1)
        "Sink"
        annotation (Placement(transformation(extent={{100,-30},{80,-10}})));
      BoilerPolynomialFurnaceHeatBalance boi(
        m1_flow_nominal=7,
        m2_flow_nominal=m_flow_nominal,
        show_T=true,
        p_start=900000,
        redeclare package MediumAir = MediumAir,
        redeclare package MediumWat = MediumWat,
        redeclare package MediumSte = MediumSte,
        m_flow_nominal=10,
        dp_nominal=600000,
        Q_flow_nominal=Q_flow_nominal,
        fue=Buildings.Fluid.Data.Fuels.Generic(
                h=47522970,
                d=0.8,
                mCO2=2.25),
        UA=0.01*Q_flow_nominal/100,
        V=2*8.76,
        V_com=20,
        FA_ratio=1.15,
        T_exh_nominal(displayUnit="K") = 439)
        annotation (Placement(transformation(extent={{8,10},{28,30}})));
      Modelica.Blocks.Sources.Ramp  ramp(height=1, duration=3600)
        annotation (Placement(transformation(extent={{-120,70},{-100,90}})));
      Modelica.Blocks.Sources.RealExpression QFue(y=boi.QFue_flow)
        "Fuel heat flow rate"
        annotation (Placement(transformation(extent={{-190,-126},{-170,-106}})));
      Modelica.Blocks.Sources.RealExpression QWat(y=boi.QWat_flow)
        "Water heat flow rate"
        annotation (Placement(transformation(extent={{-190,-176},{-170,-156}})));
      Modelica.Blocks.Sources.RealExpression QLoss(y=boi.heatPort.Q_flow)
        "Loss from boiler casing"
        annotation (Placement(transformation(extent={{-82,-178},{-62,-158}})));
      Modelica.Blocks.Sources.RealExpression QFlue(y=boi.qExhFlo.y)
        "Heat losses into flue gas"
        annotation (Placement(transformation(extent={{-82,-128},{-62,-108}})));
      Modelica.Blocks.Sources.RealExpression Qcomb(y=((boi.port_b1.h_outflow - boi.port_a1.h_outflow)
            *boi.port_a1.m_flow))
        "Q calcuated based on the mflow and enthalpy of the combustion side ports"
        annotation (Placement(transformation(extent={{46,-128},{66,-108}})));
      Modelica.Blocks.Sources.RealExpression Qevap(y=((boi.port_b2.h_outflow - boi.port_a2.h_outflow)
            *boi.port_b2.m_flow))
        "Q calcuated based on the mflow and enthalpy of the combustion side ports"
        annotation (Placement(transformation(extent={{46,-178},{66,-158}})));
      Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-152,-126},{-132,-106}})));
      Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-152,-176},{-132,-156}})));
      Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-40,-128},{-20,-108}})));
      Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-40,-178},{-20,-158}})));
      Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{94,-128},{114,-108}})));
      Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{94,-178},{114,-158}})));
      Modelica.Blocks.Math.Division perQWat
        annotation (Placement(transformation(extent={{-116,-150},{-96,-170}})));
      Modelica.Blocks.Math.Division perQLoss
        annotation (Placement(transformation(extent={{-6,-152},{14,-172}})));
      Modelica.Blocks.Math.Division perQFlue
        annotation (Placement(transformation(extent={{-4,-134},{16,-114}})));
      Modelica.Blocks.Math.Division perQcomb
        annotation (Placement(transformation(extent={{130,-134},{150,-114}})));
      Modelica.Blocks.Math.Division perQeva
        annotation (Placement(transformation(extent={{130,-152},{150,-172}})));
      Modelica.Blocks.Math.Add diffQevap(k1=-1)
        annotation (Placement(transformation(extent={{12,-206},{-8,-226}})));
      Modelica.Blocks.Math.Add diffQcomb(k1=-1)
        annotation (Placement(transformation(extent={{70,-206},{50,-226}})));
    equation
      connect(sou.ports[1],fwPum. port_a)
        annotation (Line(points={{-60,-20},{-40,-20}}, color={0,127,255},
          thickness=0.5));
      connect(mFloFW.y,fwPum. m_flow_in)
        annotation (Line(points={{-99,0},{-30,0},{-30,-8}},      color={0,0,127},
          pattern=LinePattern.Dash));
      connect(Tfw_in.y, sou.T_in) annotation (Line(points={{-97,-40},{-88,-40},
              {-88,-16},{-82,-16}},
                               color={0,0,127},
          pattern=LinePattern.Dash));
      connect(Tair_in.y, pro.T_in) annotation (Line(
          points={{-99,44},{-82,44}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(fwPum.port_b, boi.port_a2) annotation (Line(points={{-20,-20},{2,
              -20},{2,16},{8,16}},
                              color={0,127,255},
          thickness=0.5));
      connect(boi.port_b2, sin.ports[1]) annotation (Line(points={{28,16},{60,
              16},{60,-20},{80,-20}},
                              color={238,46,47},
          thickness=0.5));
      connect(boi.port_b1, exh.ports[1]) annotation (Line(
          points={{28,24},{60,24},{60,40},{80,40}},
          color={0,127,255},
          thickness=0.5));
      connect(pro.ports[1], boi.port_a1) annotation (Line(
          points={{-60,40},{-2,40},{-2,24},{8,24}},
          color={0,127,255},
          thickness=0.5));
      connect(TAmb.port, boi.heatPort) annotation (Line(points={{6,94},{18,94},
              {18,30}}, color={191,0,0}));
      connect(ramp.y, boi.y) annotation (Line(
          points={{-99,80},{-20,80},{-20,42},{0,42},{0,28},{6,28}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(QFue.y,iQfue. u)
        annotation (Line(points={{-169,-116},{-154,-116}},
                                                         color={0,0,127}));
      connect(QWat.y,iQWat. u)
        annotation (Line(points={{-169,-166},{-154,-166}}, color={0,0,127}));
      connect(QFlue.y,iQFlue. u)
        annotation (Line(points={{-61,-118},{-42,-118}},
                                                       color={0,0,127}));
      connect(QLoss.y,iQLoss. u)
        annotation (Line(points={{-61,-168},{-42,-168}}, color={0,0,127}));
      connect(iQComb.u,Qcomb. y)
        annotation (Line(points={{92,-118},{67,-118}},color={0,0,127}));
      connect(iQevap.u,Qevap. y)
        annotation (Line(points={{92,-168},{67,-168}},  color={0,0,127}));
      connect(iQWat.y,perQWat. u1)
        annotation (Line(points={{-131,-166},{-118,-166}},color={0,0,127}));
      connect(iQfue.y,perQWat. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-154},{-118,-154}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQLoss.y,perQLoss. u1)
        annotation (Line(points={{-19,-168},{-8,-168}},
                                                      color={0,0,127}));
      connect(iQfue.y,perQLoss. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-156},{-8,-156}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQFlue.y,perQFlue. u1)
        annotation (Line(points={{-19,-118},{-6,-118}},
                                                    color={0,0,127}));
      connect(iQfue.y,perQFlue. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-130},{-6,-130}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQComb.y,perQcomb. u1)
        annotation (Line(points={{115,-118},{128,-118}},
                                                       color={0,0,127}));
      connect(iQfue.y,perQcomb. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-130},{128,-130}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQevap.y,perQeva. u1)
        annotation (Line(points={{115,-168},{128,-168}}, color={0,0,127}));
      connect(iQfue.y,perQeva. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-156},{128,-156}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(perQcomb.y,diffQcomb. u2) annotation (Line(points={{151,-124},{158,-124},
              {158,-210},{72,-210}}, color={0,0,127}));
      connect(perQeva.y,diffQevap. u2) annotation (Line(points={{151,-162},{150,-162},
              {150,-202},{14,-202},{14,-210}}, color={0,0,127}));
      connect(perQWat.y,diffQevap. u1) annotation (Line(points={{-95,-160},{-88,-160},
              {-88,-236},{26,-236},{26,-222},{14,-222}}, color={244,125,35}));
      connect(perQFlue.y,diffQcomb. u1) annotation (Line(points={{17,-124},{40,-124},
              {40,-236},{82,-236},{82,-222},{72,-222}},  color={244,125,35}));
      annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Components/HeatBalanceBoiler.mos"
            "Simulate and plot"),
        experiment(Tolerance=1e-6, StopTime=3600),Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
                {160,100}})),                                        Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}}),
            graphics={
              Rectangle(
              extent={{-200,-100},{160,-240}},
              lineColor={28,108,200},
              fillColor={174,179,179},
              fillPattern=FillPattern.Solid), Text(
              extent={{-198,-82},{-92,-102}},
              textColor={28,108,200},
              textString="Steady state verification (Heat balance)")}));
    end HeatBalanceBoiler;

    model HeatBalanceBoilerCalibration
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

      Fluid.Sources.Boundary_pT           sou(
        redeclare package Medium = MediumWat,
        p=100000,
        use_T_in=true,
        T=303.15,
        nPorts=1)
        "Source"
        annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
      Fluid.Movers.FlowControlled_m_flow fwPum(
        redeclare package Medium = MediumWat,
        m_flow_nominal=m_flow_nominal,
        addPowerToMedium=false,
        nominalValuesDefineDefaultPressureCurve=true,
        dp_nominal=dp_nominal) "Feed water pump"
        annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
      Modelica.Blocks.Sources.CombiTimeTable
                                   mFloFW
        annotation (Placement(transformation(extent={{-120,-10},{-100,10}})));
      Modelica.Blocks.Sources.CombiTimeTable
                                   Tfw_in
        annotation (Placement(transformation(extent={{-120,-50},{-100,-30}})));
      Modelica.Blocks.Sources.CombiTimeTable
                                       Tair_in(offset=273)
        annotation (Placement(transformation(extent={{-120,34},{-100,54}})));
      Fluid.Sources.Boundary_pT pro(
        redeclare package Medium = MediumAir,
        use_T_in=true,
        T=573.15,
        nPorts=1) "Source"
        annotation (Placement(transformation(extent={{-80,30},{-60,50}})));
      HeatTransfer.Sources.FixedTemperature           TAmb(T(displayUnit="K") = 284)
        "Ambient temperature in boiler room"
        annotation (Placement(transformation(extent={{-14,84},{6,104}})));
      Fluid.Sources.Boundary_pT exh(
        redeclare package Medium = MediumAir,
        p=100000,
        T(displayUnit="K") = 560,
        nPorts=1) "Source"
        annotation (Placement(transformation(extent={{120,14},{100,34}})));

      Fluid.Sources.Boundary_pT           sin(
        redeclare package Medium = MediumSte,
        p(displayUnit="bar") = 900000,
        T=453.15,
        nPorts=1)
        "Sink"
        annotation (Placement(transformation(extent={{120,-30},{100,-10}})));
      BoilerPolynomialFurnaceHeatBalance boi(
        m1_flow_nominal=7,
        m2_flow_nominal=m_flow_nominal,
        show_T=true,
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        p_start=600000,
        redeclare package MediumAir = MediumAir,
        redeclare package MediumWat = MediumWat,
        redeclare package MediumSte = MediumSte,
        mDry=1.5E-3*Q_flow_nominal,
        m_flow_nominal=10,
        dp_nominal=600000,
        Q_flow_nominal=Q_flow_nominal,
        fue=Buildings.Fluid.Data.Fuels.Generic(
                h=47522970,
                d=0.8,
                mCO2=2.25),
        UA=0.01*Q_flow_nominal/100,
        V=2*8.76,
        V_com=8,
        FA_ratio=1.15,
        T_exh_nominal(displayUnit="K") = 420)
        annotation (Placement(transformation(extent={{8,10},{28,30}})));
      Modelica.Blocks.Sources.CombiTimeTable pLR(
        tableOnFile=true,
        fileName=
            "Modelica.Utilities.Files.loadResource(modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerInputs.mos)",
        columns=4:5)
        annotation (Placement(transformation(extent={{-120,70},{-100,90}})));

      Modelica.Blocks.Sources.RealExpression QFue(y=boi.QFue_flow)
        "Fuel heat flow rate"
        annotation (Placement(transformation(extent={{-190,-126},{-170,-106}})));
      Modelica.Blocks.Sources.RealExpression QWat(y=boi.QWat_flow)
        "Water heat flow rate"
        annotation (Placement(transformation(extent={{-190,-176},{-170,-156}})));
      Modelica.Blocks.Sources.RealExpression QLoss(y=boi.heatPort.Q_flow)
        "Loss from boiler casing"
        annotation (Placement(transformation(extent={{-82,-178},{-62,-158}})));
      Modelica.Blocks.Sources.RealExpression QFlue(y=boi.qExhFlo.y)
        "Heat losses into flue gas"
        annotation (Placement(transformation(extent={{-82,-128},{-62,-108}})));
      Modelica.Blocks.Sources.RealExpression Qcomb(y=((boi.port_b1.h_outflow - boi.port_a1.h_outflow)
            *boi.port_a1.m_flow))
        "Q calcuated based on the mflow and enthalpy of the combustion side ports"
        annotation (Placement(transformation(extent={{46,-128},{66,-108}})));
      Modelica.Blocks.Sources.RealExpression Qevap(y=((boi.port_b2.h_outflow - boi.port_a2.h_outflow)
            *boi.port_b2.m_flow))
        "Q calcuated based on the mflow and enthalpy of the combustion side ports"
        annotation (Placement(transformation(extent={{46,-178},{66,-158}})));
      Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-152,-126},{-132,-106}})));
      Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-152,-176},{-132,-156}})));
      Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-40,-128},{-20,-108}})));
      Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{-40,-178},{-20,-158}})));
      Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{94,-128},{114,-108}})));
      Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
        annotation (Placement(transformation(extent={{94,-178},{114,-158}})));
      Modelica.Blocks.Math.Division perQWat
        annotation (Placement(transformation(extent={{-116,-150},{-96,-170}})));
      Modelica.Blocks.Math.Division perQLoss
        annotation (Placement(transformation(extent={{-6,-152},{14,-172}})));
      Modelica.Blocks.Math.Division perQFlue
        annotation (Placement(transformation(extent={{-4,-134},{16,-114}})));
      Modelica.Blocks.Math.Division perQcomb
        annotation (Placement(transformation(extent={{130,-134},{150,-114}})));
      Modelica.Blocks.Math.Division perQeva
        annotation (Placement(transformation(extent={{130,-152},{150,-172}})));
      Modelica.Blocks.Math.Add diffQevap(k1=-1)
        annotation (Placement(transformation(extent={{12,-206},{-8,-226}})));
      Modelica.Blocks.Math.Add diffQcomb(k1=-1)
        annotation (Placement(transformation(extent={{70,-206},{50,-226}})));
      Fluid.Sensors.TemperatureTwoPort
                                     senTem(
        redeclare package Medium = MediumAir,
        m_flow_nominal=m_flow_nominal,
        tau=30,
        T_start(displayUnit="K"))
        annotation (Placement(transformation(extent={{-10,10},{10,-10}},
            rotation=180,
            origin={70,24})));
      Fluid.Sensors.TemperatureTwoPort
                                     senTem1(
        redeclare package Medium = MediumAir,
        m_flow_nominal=m_flow_nominal,
        tau=30,
        T_start(displayUnit="K"))
        annotation (Placement(transformation(extent={{-10,10},{10,-10}},
            rotation=180,
            origin={-40,40})));
      Fluid.Sensors.TemperatureTwoPort
                                     senTem2(
        redeclare package Medium = MediumSte,
        m_flow_nominal=m_flow_nominal,
        tau=30,
        T_start(displayUnit="K"))
        annotation (Placement(transformation(extent={{-10,10},{10,-10}},
            rotation=180,
            origin={70,-20})));
    equation
      connect(sou.ports[1],fwPum. port_a)
        annotation (Line(points={{-60,-20},{-40,-20}}, color={0,127,255},
          thickness=0.5));
      connect(mFloFW.y,fwPum. m_flow_in)
        annotation (Line(points={{-99,0},{-30,0},{-30,-8}},      color={0,0,127},
          pattern=LinePattern.Dash));
      connect(Tfw_in.y, sou.T_in) annotation (Line(points={{-99,-40},{-88,-40},
              {-88,-16},{-82,-16}},
                               color={0,0,127},
          pattern=LinePattern.Dash));
      connect(Tair_in.y, pro.T_in) annotation (Line(
          points={{-99,44},{-82,44}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(fwPum.port_b, boi.port_a2) annotation (Line(points={{-20,-20},{0,
              -20},{0,16},{8,16}},
                              color={0,127,255},
          thickness=0.5));
      connect(TAmb.port, boi.heatPort) annotation (Line(points={{6,94},{18,94},
              {18,30}}, color={191,0,0}));
      connect(pLR.y, boi.y) annotation (Line(
          points={{-99,80},{-20,80},{-20,42},{0,42},{0,28},{6,28}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(QFue.y,iQfue. u)
        annotation (Line(points={{-169,-116},{-154,-116}},
                                                         color={0,0,127}));
      connect(QWat.y,iQWat. u)
        annotation (Line(points={{-169,-166},{-154,-166}}, color={0,0,127}));
      connect(QFlue.y,iQFlue. u)
        annotation (Line(points={{-61,-118},{-42,-118}},
                                                       color={0,0,127}));
      connect(QLoss.y,iQLoss. u)
        annotation (Line(points={{-61,-168},{-42,-168}}, color={0,0,127}));
      connect(iQComb.u,Qcomb. y)
        annotation (Line(points={{92,-118},{67,-118}},color={0,0,127}));
      connect(iQevap.u,Qevap. y)
        annotation (Line(points={{92,-168},{67,-168}},  color={0,0,127}));
      connect(iQWat.y,perQWat. u1)
        annotation (Line(points={{-131,-166},{-118,-166}},color={0,0,127}));
      connect(iQfue.y,perQWat. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-154},{-118,-154}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQLoss.y,perQLoss. u1)
        annotation (Line(points={{-19,-168},{-8,-168}},
                                                      color={0,0,127}));
      connect(iQfue.y,perQLoss. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-156},{-8,-156}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQFlue.y,perQFlue. u1)
        annotation (Line(points={{-19,-118},{-6,-118}},
                                                    color={0,0,127}));
      connect(iQfue.y,perQFlue. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-130},{-6,-130}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQComb.y,perQcomb. u1)
        annotation (Line(points={{115,-118},{128,-118}},
                                                       color={0,0,127}));
      connect(iQfue.y,perQcomb. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-130},{128,-130}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(iQevap.y,perQeva. u1)
        annotation (Line(points={{115,-168},{128,-168}}, color={0,0,127}));
      connect(iQfue.y,perQeva. u2) annotation (Line(
          points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-156},{128,-156}},
          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(perQcomb.y,diffQcomb. u2) annotation (Line(points={{151,-124},{158,-124},
              {158,-210},{72,-210}}, color={0,0,127}));
      connect(perQeva.y,diffQevap. u2) annotation (Line(points={{151,-162},{150,-162},
              {150,-202},{14,-202},{14,-210}}, color={0,0,127}));
      connect(perQWat.y,diffQevap. u1) annotation (Line(points={{-95,-160},{-88,-160},
              {-88,-236},{26,-236},{26,-222},{14,-222}}, color={244,125,35}));
      connect(perQFlue.y,diffQcomb. u1) annotation (Line(points={{17,-124},{40,-124},
              {40,-236},{82,-236},{82,-222},{72,-222}},  color={244,125,35}));
      connect(boi.port_b1, senTem.port_b) annotation (Line(
          points={{28,24},{60,24}},
          color={0,0,0},
          thickness=0.5));
      connect(senTem.port_a, exh.ports[1]) annotation (Line(
          points={{80,24},{100,24}},
          color={0,0,0},
          thickness=0.5));
      connect(pro.ports[1], senTem1.port_b) annotation (Line(
          points={{-60,40},{-50,40}},
          color={0,0,0},
          thickness=0.5));
      connect(senTem1.port_a, boi.port_a1) annotation (Line(
          points={{-30,40},{-2,40},{-2,24},{8,24}},
          color={0,0,0},
          thickness=0.5));
      connect(boi.port_b2, senTem2.port_b) annotation (Line(
          points={{28,16},{40,16},{40,-20},{60,-20}},
          color={238,46,47},
          thickness=0.5));
      connect(senTem2.port_a, sin.ports[1]) annotation (Line(
          points={{80,-20},{100,-20}},
          color={238,46,47},
          thickness=0.5));
      annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Components/HeatBalanceBoiler.mos"
            "Simulate and plot"),
        experiment(Tolerance=1e-6, StopTime=3600),Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
                {160,100}})),                                        Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}}),
            graphics={
              Rectangle(
              extent={{-200,-100},{160,-240}},
              lineColor={28,108,200},
              fillColor={174,179,179},
              fillPattern=FillPattern.Solid), Text(
              extent={{-198,-82},{-92,-102}},
              textColor={28,108,200},
              textString="Steady state verification (Heat balance)"),
              Rectangle(
              extent={{160,100},{300,-240}},
              lineColor={28,108,200},
              fillColor={174,179,179},
              fillPattern=FillPattern.Solid)}));
    end HeatBalanceBoilerCalibration;

    package Calibration

      model HeatBalanceBoilerTrain
          extends Modelica.Icons.Example;

        // Medium declarations
        package MediumWat =
            Buildings.Media.Specialized.Water.TemperatureDependentDensity
          "Water medium - port_a (inlet)";
        package MediumSte = Buildings.Media.Steam (
          p_default=900000,
          T_default=273.15+180,
          h_default=2777170)
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

             parameter String Inputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTrain.mos");
           parameter String Outputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTrainOut.mos");

        Fluid.Sources.Boundary_pT           sou(
          redeclare package Medium = MediumWat,
          p=100000,
          use_T_in=true,
          T=303.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
        Fluid.Movers.FlowControlled_m_flow fwPum(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m_flow_nominal,
          addPowerToMedium=false,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal=dp_nominal) "Feed water pump"
          annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
        Fluid.Sources.Boundary_pT pro(
          redeclare package Medium = MediumAir,
          use_T_in=true,
          T=573.15,
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{-80,30},{-60,50}})));
        HeatTransfer.Sources.FixedTemperature           TAmb(T(displayUnit="K") = 305)
          "Ambient temperature in boiler room"
          annotation (Placement(transformation(extent={{-14,84},{6,104}})));
        Fluid.Sources.Boundary_pT exh(
          redeclare package Medium = MediumAir,
          T(displayUnit="K"),
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{150,14},{130,34}})));

        Fluid.Sources.Boundary_pT           sin(
          redeclare package Medium = MediumSte,
          p(displayUnit="bar") = 895000,
          T=453.15,
          nPorts=1)
          "Sink"
          annotation (Placement(transformation(extent={{150,-30},{130,-10}})));
        BoilerPolynomialFurnaceHeatBalance boi(
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
          annotation (Placement(transformation(extent={{8,10},{28,30}})));

        Modelica.Blocks.Sources.RealExpression QFue(y=boi.QFue_flow)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{-190,-126},{-170,-106}})));
        Modelica.Blocks.Sources.RealExpression QWat(y=boi.QWat_flow)
          "Water heat flow rate"
          annotation (Placement(transformation(extent={{-190,-176},{-170,-156}})));
        Modelica.Blocks.Sources.RealExpression QLoss(y=boi.heatPort.Q_flow)
          "Loss from boiler casing"
          annotation (Placement(transformation(extent={{-82,-178},{-62,-158}})));
        Modelica.Blocks.Sources.RealExpression QFlue(y=boi.qExhFlo.y)
          "Heat losses into flue gas"
          annotation (Placement(transformation(extent={{-82,-128},{-62,-108}})));
        Modelica.Blocks.Sources.RealExpression Qcomb(y=((boi.port_b1.h_outflow - boi.port_a1.h_outflow)
              *boi.port_a1.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{46,-128},{66,-108}})));
        Modelica.Blocks.Sources.RealExpression Qevap(y=((boi.port_b2.h_outflow - boi.port_a2.h_outflow)
              *boi.port_b2.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{46,-178},{66,-158}})));
        Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-152,-126},{-132,-106}})));
        Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-152,-176},{-132,-156}})));
        Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-40,-128},{-20,-108}})));
        Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-40,-178},{-20,-158}})));
        Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{94,-128},{114,-108}})));
        Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{94,-178},{114,-158}})));
        Modelica.Blocks.Math.Division perQWat
          annotation (Placement(transformation(extent={{-116,-150},{-96,-170}})));
        Modelica.Blocks.Math.Division perQLoss
          annotation (Placement(transformation(extent={{-6,-152},{14,-172}})));
        Modelica.Blocks.Math.Division perQFlue
          annotation (Placement(transformation(extent={{-4,-134},{16,-114}})));
        Modelica.Blocks.Math.Division perQcomb
          annotation (Placement(transformation(extent={{130,-134},{150,-114}})));
        Modelica.Blocks.Math.Division perQeva
          annotation (Placement(transformation(extent={{130,-152},{150,-172}})));
        Modelica.Blocks.Math.Add diffQevap(k1=-1)
          annotation (Placement(transformation(extent={{12,-206},{-8,-226}})));
        Modelica.Blocks.Math.Add diffQcomb(k1=-1)
          annotation (Placement(transformation(extent={{70,-206},{50,-226}})));
        Fluid.Sensors.TemperatureTwoPort
                                       senTem(
          redeclare package Medium = MediumAir,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K"))
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=180,
              origin={70,24})));
        Fluid.Sensors.TemperatureTwoPort
                                       senTem1(
          redeclare package Medium = MediumAir,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K"))
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=180,
              origin={-40,40})));
        Fluid.Sensors.TemperatureTwoPort
                                       senTem2(
          redeclare package Medium = MediumSte,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K"))
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=180,
              origin={70,-20})));
        Modelica.Blocks.Sources.RealExpression Qboi_S(y=QWat.y - QLoss.y)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{176,28},{196,48}})));
        Modelica.Blocks.Sources.RealExpression mFueFlo_S(y=boi.mFue_flow)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,8},{196,28}})));
        Modelica.Blocks.Sources.RealExpression T_exh_S(y=senTem.T)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,-12},{196,8}})));
        Modelica.Blocks.Sources.RealExpression T_Ste_S(y=senTem2.T)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,-32},{196,-12}})));
        Modelica.Blocks.Sources.CombiTimeTable pLR(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={5},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-120,70},{-100,90}})));
        Modelica.Blocks.Sources.CombiTimeTable TAirIn(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={2},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-180,34},{-160,54}})));
        Modelica.Blocks.Sources.CombiTimeTable mFloFw(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={3},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-180,0},{-160,20}})));
        Modelica.Blocks.Sources.CombiTimeTable TFw(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={4},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-178,-42},{-158,-22}})));
        Modelica.Blocks.Sources.CombiTimeTable MeaData(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Outputs),
          verboseRead=true,
          columns=2:5,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{176,62},{196,82}})));
        Fluid.FixedResistances.PressureDrop res(
          redeclare package Medium = MediumAir,
          m_flow_nominal=10,
          dp_nominal=50000)
          annotation (Placement(transformation(extent={{100,14},{120,34}})));
        Fluid.FixedResistances.PressureDrop res1(
          redeclare package Medium = MediumSte,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=10000)
          annotation (Placement(transformation(extent={{100,-30},{120,-10}})));
        Modelica.Blocks.Sources.CombiTimeTable stePre(
          tableOnFile=false,
          table=[8683200,892316; 8684100,891685; 8685000,890902; 8685900,891196; 8686800,
              895785; 8687700,898587; 8688600,900129; 8689500,896132; 8690400,891488;
              8691300,885147; 8692200,889115; 8693100,896804; 8694000,893937; 8694900,
              894953; 8695800,892835; 8696700,891911; 8697600,893801; 8698500,897013;
              8699400,890919; 8700300,894550; 8701200,894498; 8702100,891167; 8703000,
              887755; 8703900,891254; 8704800,891116; 8705700,895712; 8706600,891681;
              8707500,890279; 8708400,893290; 8709300,895085; 8710200,893762; 8711100,
              893850; 8712000,891206; 8712900,889010; 8713800,891960; 8714700,893310;
              8715600,894229; 8716500,894933; 8717400,892706; 8718300,890688; 8719200,
              892198; 8720100,888249; 8721000,890422; 8721900,895005; 8722800,894123;
              8723700,893963; 8724600,892238; 8725500,891018; 8726400,893750; 8727300,
              892906; 8728200,891890; 8729100,892280; 8730000,890381; 8730900,889701;
              8731800,892948; 8732700,890220; 8733600,892233; 8734500,891652; 8735400,
              888630; 8736300,890460; 8737200,890518; 8738100,890228; 8739000,895814;
              8739900,900284; 8740800,900628; 8741700,898741; 8742600,898387; 8743500,
              896780; 8744400,895804; 8745300,896948; 8746200,896786; 8747100,897274;
              8748000,892660; 8748900,897156; 8749800,896337; 8750700,897688; 8751600,
              892595; 8752500,897295; 8753400,898425; 8754300,894919; 8755200,899118;
              8756100,896767; 8757000,897065; 8757900,894714; 8758800,889521; 8759700,
              901560; 8760600,900437; 8761500,901380; 8762400,898418; 8763300,899716;
              8764200,897777; 8765100,900712; 8766000,898616; 8766900,898191; 8767800,
              897273; 8768700,900192; 8769600,902553; 8770500,901377; 8771400,901969;
              8772300,899335; 8773200,900886; 8774100,899350; 8775000,900286; 8775900,
              898729; 8776800,902840; 8777700,899746; 8778600,900780; 8779500,901832;
              8780400,901287; 8781300,901840; 8782200,901655; 8783100,900987; 8784000,
              899158; 8784900,899444; 8785800,897361; 8786700,897411; 8787600,896690;
              8788500,894345; 8789400,897702; 8790300,899003; 8791200,899007; 8792100,
              897467; 8793000,898909; 8793900,896982; 8794800,899785; 8795700,898488;
              8796600,897608; 8797500,901257; 8798400,898667; 8799300,897137; 8800200,
              899305; 8801100,899503; 8802000,899118; 8802900,900084; 8803800,898564;
              8804700,899755; 8805600,902665; 8806500,898058; 8807400,901453; 8808300,
              900676; 8809200,897994; 8810100,902194; 8811000,898462; 8811900,897836;
              8812800,902642; 8813700,902053; 8814600,898060; 8815500,898509; 8816400,
              901198; 8817300,902410; 8818200,900042; 8819100,901264; 8820000,901855;
              8820900,903825; 8821800,898006; 8822700,902393; 8823600,897728; 8824500,
              897881; 8825400,899065; 8826300,899409; 8827200,900042; 8828100,899106;
              8829000,894843; 8829900,897442; 8830800,891668; 8831700,899983; 8832600,
              900388; 8833500,898621; 8834400,894455; 8835300,899183; 8836200,902382;
              8837100,901714; 8838000,892222; 8838900,897753; 8839800,895468; 8840700,
              896568; 8841600,900810; 8842500,898862; 8843400,902657; 8844300,907208;
              8845200,901789; 8846100,905161; 8847000,904069; 8847900,919049; 8848800,
              933409; 8849700,933443; 8850600,937552; 8851500,917675; 8852400,914861;
              8853300,927014; 8854200,924750; 8855100,920426; 8856000,909951; 8856900,
              913132; 8857800,920944; 8858700,905232; 8859600,904585; 8860500,903512;
              8861400,903256; 8862300,899210; 8863200,899369; 8864100,897609; 8865000,
              899453; 8865900,899356; 8866800,898987; 8867700,898571; 8868600,897670;
              8869500,898591; 8870400,898527; 8871300,897947; 8872200,892801; 8873100,
              891317; 8874000,898983; 8874900,892916; 8875800,897873; 8876700,896933;
              8877600,895837; 8878500,896179; 8879400,900107; 8880300,897665; 8881200,
              898789; 8882100,906752; 8883000,898733; 8883900,896273; 8884800,903821;
              8885700,896866; 8886600,898431; 8887500,901278; 8888400,882846; 8889300,
              895956; 8890200,902357; 8891100,901307; 8892000,899826; 8892900,900040;
              8893800,896459; 8894700,897382; 8895600,899511; 8896500,896183; 8897400,
              903493; 8898300,899879; 8899200,899473; 8900100,902817; 8901000,897306;
              8901900,901992; 8902800,900536; 8903700,903393; 8904600,901198; 8905500,
              902132; 8906400,899723; 8907300,897768; 8908200,898384; 8909100,894115;
              8910000,894856; 8910900,892514; 8911800,901170; 8912700,896023; 8913600,
              901879; 8914500,897391; 8915400,899412; 8916300,904120; 8917200,899862;
              8918100,894087; 8919000,895448; 8919900,895685; 8920800,909290; 8921700,
              905152; 8922600,895080; 8923500,900401; 8924400,902049; 8925300,897232;
              8926200,895718; 8927100,900294; 8928000,901106; 8928900,901236; 8929800,
              902719; 8930700,900085; 8931600,902612; 8932500,901231; 8933400,902378;
              8934300,884259; 8935200,888231; 8936100,899815; 8937000,893617; 8937900,
              896661; 8938800,902966; 8939700,908549; 8940600,901597; 8941500,894364;
              8942400,895806; 8943300,901993; 8944200,893605; 8945100,887119; 8946000,
              894041; 8946900,898020; 8947800,898204; 8948700,901734; 8949600,906864;
              8950500,896506; 8951400,895502; 8952300,890127; 8953200,902039; 8954100,
              902657; 8955000,881041; 8955900,881701; 8956800,907420; 8957700,901602;
              8958600,900108; 8959500,895231; 8960400,896903; 8961300,899413; 8962200,
              895501; 8963100,895750; 8964000,898890; 8964900,897228; 8965800,896884;
              8966700,896233; 8967600,897827; 8968500,896559; 8969400,896825; 8970300,
              896490; 8971200,898713; 8972100,896536; 8973000,895528; 8973900,896828;
              8974800,897874; 8975700,897805; 8976600,900293; 8977500,895008; 8978400,
              894395; 8979300,897377; 8980200,896797; 8981100,901117; 8982000,901607;
              8982900,896163; 8983800,890473; 8984700,895937; 8985600,900544; 8986500,
              899973; 8987400,901876; 8988300,899390; 8989200,896199; 8990100,898284;
              8991000,899128; 8991900,901099; 8992800,898263; 8993700,896992; 8994600,
              900369; 8995500,899205; 8996400,901288; 8997300,903209; 8998200,899505;
              8999100,898207; 9000000,889801; 9000900,893753; 9001800,900226; 9002700,
              899721; 9003600,902705; 9004500,904361; 9005400,902308; 9006300,889626;
              9007200,885060; 9008100,894836; 9009000,903038; 9009900,894932; 9010800,
              893072; 9011700,906769; 9012600,904800; 9013500,901890; 9014400,900425;
              9015300,907985; 9016200,911117; 9017100,913138; 9018000,921595; 9018900,
              938608; 9019800,952207; 9020700,923140; 9021600,900445; 9022500,905692;
              9023400,908128; 9024300,889486; 9025200,897555; 9026100,896083; 9027000,
              888324; 9027900,889318; 9028800,889185; 9029700,898645; 9030600,899778;
              9031500,897199; 9032400,898766; 9033300,898734; 9034200,896020; 9035100,
              891446; 9036000,899859; 9036900,899089; 9037800,897368; 9038700,898364;
              9039600,896495; 9040500,894885; 9041400,894342; 9042300,899818; 9043200,
              899519; 9044100,896962; 9045000,895953; 9045900,894465; 9046800,895224;
              9047700,894825; 9048600,895060; 9049500,897563; 9050400,894841; 9051300,
              894129; 9052200,897354; 9053100,897693; 9054000,899312; 9054900,895058;
              9055800,894160; 9056700,896299; 9057600,896457; 9058500,897513; 9059400,
              898402; 9060300,897171; 9061200,895457; 9062100,897859; 9063000,897867;
              9063900,897006; 9064800,899385; 9065700,895948; 9066600,894770; 9067500,
              897424; 9068400,894896; 9069300,892490; 9070200,890507; 9071100,894937;
              9072000,897981; 9072900,898133; 9073800,901646; 9074700,901137; 9075600,
              896753; 9076500,890813; 9077400,896331; 9078300,896888; 9079200,904427;
              9080100,901888; 9081000,893080; 9081900,891112; 9082800,895399; 9083700,
              897208; 9084600,900491; 9085500,894977; 9086400,893689; 9087300,893570;
              9088200,890527; 9089100,894471; 9090000,891748; 9090900,894354; 9091800,
              896214; 9092700,894799; 9093600,890308; 9094500,893537; 9095400,893195;
              9096300,890145; 9097200,891324; 9098100,886895; 9099000,888348; 9099900,
              902638; 9100800,896894; 9101700,894346; 9102600,893451; 9103500,893956;
              9104400,904706; 9105300,900887; 9106200,895223; 9107100,892098; 9108000,
              898905; 9108900,896420; 9109800,896890; 9110700,900077; 9111600,901466;
              9112500,898862; 9113400,900364; 9114300,899451],
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{176,-66},{196,-46}})));
        Modelica.Blocks.Sources.RealExpression P_Ste_S(y=boi.port_b2.p)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{178,-92},{198,-72}})));
        Modelica.Blocks.Math.Gain             iQfue1(k=0.8)
          annotation (Placement(transformation(extent={{-140,0},{-120,20}})));
      equation
        connect(sou.ports[1],fwPum. port_a)
          annotation (Line(points={{-60,-20},{-40,-20}}, color={0,127,255},
            thickness=0.5));
        connect(fwPum.port_b, boi.port_a2) annotation (Line(points={{-20,-20},{0,
                -20},{0,16},{8,16}},
                                color={0,127,255},
            thickness=0.5));
        connect(TAmb.port, boi.heatPort) annotation (Line(points={{6,94},{18,94},
                {18,30}}, color={191,0,0}));
        connect(QFue.y,iQfue. u)
          annotation (Line(points={{-169,-116},{-154,-116}},
                                                           color={0,0,127}));
        connect(QWat.y,iQWat. u)
          annotation (Line(points={{-169,-166},{-154,-166}}, color={0,0,127}));
        connect(QFlue.y,iQFlue. u)
          annotation (Line(points={{-61,-118},{-42,-118}},
                                                         color={0,0,127}));
        connect(QLoss.y,iQLoss. u)
          annotation (Line(points={{-61,-168},{-42,-168}}, color={0,0,127}));
        connect(iQComb.u,Qcomb. y)
          annotation (Line(points={{92,-118},{67,-118}},color={0,0,127}));
        connect(iQevap.u,Qevap. y)
          annotation (Line(points={{92,-168},{67,-168}},  color={0,0,127}));
        connect(iQWat.y,perQWat. u1)
          annotation (Line(points={{-131,-166},{-118,-166}},color={0,0,127}));
        connect(iQfue.y,perQWat. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-154},{-118,-154}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQLoss.y,perQLoss. u1)
          annotation (Line(points={{-19,-168},{-8,-168}},
                                                        color={0,0,127}));
        connect(iQfue.y,perQLoss. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-156},{-8,-156}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQFlue.y,perQFlue. u1)
          annotation (Line(points={{-19,-118},{-6,-118}},
                                                      color={0,0,127}));
        connect(iQfue.y,perQFlue. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-130},{-6,-130}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQComb.y,perQcomb. u1)
          annotation (Line(points={{115,-118},{128,-118}},
                                                         color={0,0,127}));
        connect(iQfue.y,perQcomb. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-130},{128,-130}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQevap.y,perQeva. u1)
          annotation (Line(points={{115,-168},{128,-168}}, color={0,0,127}));
        connect(iQfue.y,perQeva. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-156},{128,-156}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(perQcomb.y,diffQcomb. u2) annotation (Line(points={{151,-124},{158,-124},
                {158,-210},{72,-210}}, color={0,0,127}));
        connect(perQeva.y,diffQevap. u2) annotation (Line(points={{151,-162},{150,-162},
                {150,-202},{14,-202},{14,-210}}, color={0,0,127}));
        connect(perQWat.y,diffQevap. u1) annotation (Line(points={{-95,-160},{-88,-160},
                {-88,-236},{26,-236},{26,-222},{14,-222}}, color={244,125,35}));
        connect(perQFlue.y,diffQcomb. u1) annotation (Line(points={{17,-124},{40,-124},
                {40,-236},{82,-236},{82,-222},{72,-222}},  color={244,125,35}));
        connect(boi.port_b1, senTem.port_b) annotation (Line(
            points={{28,24},{60,24}},
            color={0,0,0},
            thickness=0.5));
        connect(pro.ports[1], senTem1.port_b) annotation (Line(
            points={{-60,40},{-50,40}},
            color={0,0,0},
            thickness=0.5));
        connect(senTem1.port_a, boi.port_a1) annotation (Line(
            points={{-30,40},{-2,40},{-2,24},{8,24}},
            color={0,0,0},
            thickness=0.5));
        connect(boi.port_b2, senTem2.port_b) annotation (Line(
            points={{28,16},{40,16},{40,-20},{60,-20}},
            color={238,46,47},
            thickness=0.5));
        connect(pLR.y[1], boi.y) annotation (Line(
            points={{-99,80},{-20,80},{-20,52},{0,52},{0,28},{6,28}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(senTem.port_a, res.port_a)
          annotation (Line(points={{80,24},{100,24}}, color={0,127,255}));
        connect(res.port_b, exh.ports[1])
          annotation (Line(points={{120,24},{130,24}}, color={0,127,255}));
        connect(senTem2.port_a, res1.port_a)
          annotation (Line(points={{80,-20},{100,-20}}, color={0,127,255}));
        connect(res1.port_b, sin.ports[1])
          annotation (Line(points={{120,-20},{130,-20}}, color={0,127,255}));
        connect(iQfue1.y, fwPum.m_flow_in)
          annotation (Line(points={{-119,10},{-30,10},{-30,-8}}, color={0,0,127}));
        connect(mFloFw.y[1], iQfue1.u)
          annotation (Line(points={{-159,10},{-142,10}}, color={0,0,127}));
        connect(TAirIn.y[1], pro.T_in)
          annotation (Line(points={{-159,44},{-82,44}}, color={0,0,127}));
        connect(TFw.y[1], sou.T_in) annotation (Line(points={{-157,-32},{-90,
                -32},{-90,-16},{-82,-16}}, color={0,0,127}));
        annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Calibration/BoilerFurnace.mos"
              "Simulate and plot"),
          experiment(
            StartTime=8683200,
            StopTime=9114336,
            Tolerance=1e-06,
            __Dymola_Algorithm="Dassl"),            Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
                  {160,100}})),                                        Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}}),
              graphics={
                Rectangle(
                extent={{-200,-100},{160,-240}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid), Text(
                extent={{-198,-82},{-92,-102}},
                textColor={28,108,200},
                textString="Steady state verification (Heat balance)"),
                Rectangle(
                extent={{160,100},{300,-240}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid)}));
      end HeatBalanceBoilerTrain;

      model HeatBalanceBoilerTest
          extends Modelica.Icons.Example;

        // Medium declarations
        package MediumWat =
            Buildings.Media.Specialized.Water.TemperatureDependentDensity
          "Water medium - port_a (inlet)";
        package MediumSte = Buildings.Media.Steam (
          p_default=900000,
          T_default=273.15+180,
          h_default=2777170)
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

             parameter String Inputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTestIn.mos");
           parameter String Outputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTestOut.mos");

        Fluid.Sources.Boundary_pT           sou(
          redeclare package Medium = MediumWat,
          p=100000,
          use_T_in=true,
          T=303.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
        Fluid.Movers.FlowControlled_m_flow fwPum(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m_flow_nominal,
          addPowerToMedium=false,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal=dp_nominal) "Feed water pump"
          annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
        Fluid.Sources.Boundary_pT pro(
          redeclare package Medium = MediumAir,
          use_T_in=true,
          T=573.15,
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{-80,30},{-60,50}})));
        HeatTransfer.Sources.FixedTemperature           TAmb(T(displayUnit="K") = 305)
          "Ambient temperature in boiler room"
          annotation (Placement(transformation(extent={{-14,84},{6,104}})));
        Fluid.Sources.Boundary_pT exh(
          redeclare package Medium = MediumAir,
          T(displayUnit="K"),
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{150,14},{130,34}})));

        Fluid.Sources.Boundary_pT           sin(
          redeclare package Medium = MediumSte,
          p(displayUnit="bar") = 900000,
          T=453.15,
          nPorts=1)
          "Sink"
          annotation (Placement(transformation(extent={{150,-30},{130,-10}})));
        BoilerPolynomialFurnaceHeatBalance boi(
          m1_flow_nominal=10,
          m2_flow_nominal=m_flow_nominal,
          show_T=true,
          energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
          p_start=900000,
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
          FA_ratio=1.15,
          T_exh_nominal(displayUnit="K") = 420)
          annotation (Placement(transformation(extent={{8,10},{28,30}})));

        Modelica.Blocks.Sources.RealExpression QFue(y=boi.QFue_flow)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{-190,-126},{-170,-106}})));
        Modelica.Blocks.Sources.RealExpression QWat(y=boi.QWat_flow)
          "Water heat flow rate"
          annotation (Placement(transformation(extent={{-190,-176},{-170,-156}})));
        Modelica.Blocks.Sources.RealExpression QLoss(y=boi.heatPort.Q_flow)
          "Loss from boiler casing"
          annotation (Placement(transformation(extent={{-82,-178},{-62,-158}})));
        Modelica.Blocks.Sources.RealExpression QFlue(y=boi.qExhFlo.y)
          "Heat losses into flue gas"
          annotation (Placement(transformation(extent={{-82,-128},{-62,-108}})));
        Modelica.Blocks.Sources.RealExpression Qcomb(y=((boi.port_b1.h_outflow - boi.port_a1.h_outflow)
              *boi.port_a1.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{46,-128},{66,-108}})));
        Modelica.Blocks.Sources.RealExpression Qevap(y=((boi.port_b2.h_outflow - boi.port_a2.h_outflow)
              *boi.port_b2.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{46,-178},{66,-158}})));
        Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-152,-126},{-132,-106}})));
        Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-152,-176},{-132,-156}})));
        Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-40,-128},{-20,-108}})));
        Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-40,-178},{-20,-158}})));
        Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{94,-128},{114,-108}})));
        Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{94,-178},{114,-158}})));
        Modelica.Blocks.Math.Division perQWat
          annotation (Placement(transformation(extent={{-116,-150},{-96,-170}})));
        Modelica.Blocks.Math.Division perQLoss
          annotation (Placement(transformation(extent={{-6,-152},{14,-172}})));
        Modelica.Blocks.Math.Division perQFlue
          annotation (Placement(transformation(extent={{-4,-134},{16,-114}})));
        Modelica.Blocks.Math.Division perQcomb
          annotation (Placement(transformation(extent={{130,-134},{150,-114}})));
        Modelica.Blocks.Math.Division perQeva
          annotation (Placement(transformation(extent={{130,-152},{150,-172}})));
        Modelica.Blocks.Math.Add diffQevap(k1=-1)
          annotation (Placement(transformation(extent={{12,-206},{-8,-226}})));
        Modelica.Blocks.Math.Add diffQcomb(k1=-1)
          annotation (Placement(transformation(extent={{70,-206},{50,-226}})));
        Fluid.Sensors.TemperatureTwoPort
                                       senTem(
          redeclare package Medium = MediumAir,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K"))
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=180,
              origin={70,24})));
        Fluid.Sensors.TemperatureTwoPort
                                       senTem1(
          redeclare package Medium = MediumAir,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K"))
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=180,
              origin={-40,40})));
        Fluid.Sensors.TemperatureTwoPort
                                       senTem2(
          redeclare package Medium = MediumSte,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K"))
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=180,
              origin={70,-20})));
        Modelica.Blocks.Sources.RealExpression Qboi_S(y=QWat.y - QLoss.y)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{176,28},{196,48}})));
        Modelica.Blocks.Sources.RealExpression mFueFlo_S(y=boi.mFue_flow)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,8},{196,28}})));
        Modelica.Blocks.Sources.RealExpression T_exh_S(y=senTem.T)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,-12},{196,8}})));
        Modelica.Blocks.Sources.RealExpression T_Ste_S(y=senTem2.T)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,-32},{196,-12}})));
        Modelica.Blocks.Sources.CombiTimeTable pLR(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={5},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-120,70},{-100,90}})));
        Modelica.Blocks.Sources.CombiTimeTable TAirIn(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={2},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-180,34},{-160,54}})));
        Modelica.Blocks.Sources.CombiTimeTable mFloFw(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={3},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-180,0},{-160,20}})));
        Modelica.Blocks.Sources.CombiTimeTable TFw(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={4},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-178,-42},{-158,-22}})));
        Modelica.Blocks.Sources.CombiTimeTable MeaData(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Outputs),
          verboseRead=true,
          columns=2:5,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{176,62},{196,82}})));
        Fluid.FixedResistances.PressureDrop res(
          redeclare package Medium = MediumAir,
          m_flow_nominal=10,
          dp_nominal=50000)
          annotation (Placement(transformation(extent={{100,14},{120,34}})));
        Fluid.FixedResistances.PressureDrop res1(
          redeclare package Medium = MediumSte,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=10000)
          annotation (Placement(transformation(extent={{100,-30},{120,-10}})));
        Modelica.Blocks.Sources.CombiTimeTable stePre(
          tableOnFile=false,
          table=[30499200,919414; 30500100,909991; 30501000,907472; 30501900,
              913063; 30502800,904576; 30503700,913630; 30504600,908176;
              30505500,911337; 30506400,911394; 30507300,904142; 30508200,
              913410; 30509100,910303; 30510000,908598; 30510900,910015;
              30511800,907132; 30512700,904551; 30513600,906809; 30514500,
              914276; 30515400,908386; 30516300,910500; 30517200,909198;
              30518100,909348; 30519000,909194; 30519900,907722; 30520800,
              907054; 30521700,911150; 30522600,911070; 30523500,910175;
              30524400,912450; 30525300,911907; 30526200,911925; 30527100,
              909713; 30528000,916371; 30528900,915565; 30529800,906337;
              30530700,913331; 30531600,908646; 30532500,910214; 30533400,
              908381; 30534300,913623; 30535200,908608; 30536100,910235;
              30537000,912021; 30537900,909838; 30538800,911023; 30539700,
              907088; 30540600,913307; 30541500,911445; 30542400,908662;
              30543300,911802; 30544200,910537; 30545100,911429; 30546000,
              908790; 30546900,909965; 30547800,910747; 30548700,910535;
              30549600,909052; 30550500,913429; 30551400,909303; 30552300,
              908422; 30553200,911890; 30554100,906685; 30555000,913634;
              30555900,906922; 30556800,913763; 30557700,907316; 30558600,
              912966; 30559500,906355; 30560400,911276; 30561300,908695;
              30562200,912993; 30563100,909080; 30564000,911166; 30564900,
              906257; 30565800,913348; 30566700,909264; 30567600,909067;
              30568500,911054; 30569400,907340; 30570300,912528; 30571200,
              909215; 30572100,910275; 30573000,908964; 30573900,910707;
              30574800,907988; 30575700,911548; 30576600,909278; 30577500,
              911159; 30578400,909961; 30579300,911645; 30580200,908741;
              30581100,911200; 30582000,908817; 30582900,910578; 30583800,
              903567; 30584700,911165; 30585600,913243; 30586500,904381;
              30587400,909569; 30588300,915403; 30589200,908984; 30590100,
              911976; 30591000,909001; 30591900,910139; 30592800,911785;
              30593700,908654; 30594600,904263; 30595500,913854; 30596400,
              913221; 30597300,907629; 30598200,911394; 30599100,906895;
              30600000,906380; 30600900,912758; 30601800,910516; 30602700,
              912233; 30603600,902173; 30604500,913503; 30605400,905116;
              30606300,911866; 30607200,910138; 30608100,909290; 30609000,
              910838; 30609900,910075; 30610800,910628; 30611700,907682;
              30612600,913834; 30613500,908237; 30614400,910808; 30615300,
              911302; 30616200,912769; 30617100,911523; 30618000,909228;
              30618900,911971; 30619800,910678; 30620700,911611; 30621600,
              907990; 30622500,912272; 30623400,910228; 30624300,912202;
              30625200,911967; 30626100,908153; 30627000,909324; 30627900,
              912093; 30628800,910997; 30629700,905399; 30630600,915864;
              30631500,907357; 30632400,912020; 30633300,909104; 30634200,
              909634; 30635100,910011; 30636000,912473; 30636900,910331;
              30637800,906842; 30638700,915139; 30639600,903306; 30640500,
              914123; 30641400,912408; 30642300,908970; 30643200,909106;
              30644100,912544; 30645000,907720; 30645900,911959; 30646800,
              909188; 30647700,910441; 30648600,908697; 30649500,911103;
              30650400,910377; 30651300,908953; 30652200,909650; 30653100,
              909411; 30654000,906314; 30654900,912272; 30655800,909724;
              30656700,911456; 30657600,908931; 30658500,911144; 30659400,
              905616; 30660300,915738; 30661200,907338; 30662100,911444;
              30663000,910290; 30663900,910093; 30664800,909923; 30665700,
              909880; 30666600,910289; 30667500,911092; 30668400,909928;
              30669300,909279; 30670200,909244; 30671100,911616; 30672000,
              909669; 30672900,912281; 30673800,907716; 30674700,910556;
              30675600,910154; 30676500,907298; 30677400,912371; 30678300,
              908687; 30679200,910674; 30680100,910584; 30681000,911617;
              30681900,909341; 30682800,908289; 30683700,912311; 30684600,
              908446; 30685500,911128; 30686400,906441; 30687300,912439;
              30688200,909286; 30689100,910846; 30690000,909365; 30690900,
              906423; 30691800,913213; 30692700,906253; 30693600,901091;
              30694500,915173; 30695400,914870; 30696300,903303; 30697200,
              913041; 30698100,911415; 30699000,907880; 30699900,912214;
              30700800,910461; 30701700,911762; 30702600,911261; 30703500,
              907134; 30704400,912935; 30705300,910535; 30706200,910501;
              30707100,910016; 30708000,911268; 30708900,909825; 30709800,
              913636; 30710700,911614; 30711600,908795; 30712500,898808;
              30713400,919566; 30714300,908594; 30715200,911201; 30716100,
              911736; 30717000,908918; 30717900,913010; 30718800,910472;
              30719700,907986; 30720600,907799; 30721500,911371; 30722400,
              910922; 30723300,910446; 30724200,907521; 30725100,910294;
              30726000,911425; 30726900,909628; 30727800,909927; 30728700,
              911108; 30729600,911439; 30730500,906469; 30731400,911149;
              30732300,912181; 30733200,909104; 30734100,909450; 30735000,
              911211; 30735900,908413; 30736800,909488; 30737700,914269;
              30738600,906264; 30739500,910419; 30740400,906188; 30741300,
              912745; 30742200,912658; 30743100,906504; 30744000,912897;
              30744900,909549; 30745800,908794; 30746700,912346; 30747600,
              908972; 30748500,910575; 30749400,910187; 30750300,910601;
              30751200,910459; 30752100,909811; 30753000,909424; 30753900,
              911877; 30754800,906643; 30755700,911954; 30756600,911430;
              30757500,909083; 30758400,910656; 30759300,910877; 30760200,
              907081; 30761100,911224; 30762000,912579; 30762900,910058;
              30763800,908160; 30764700,911417; 30765600,910075; 30766500,
              904861; 30767400,915895; 30768300,908200; 30769200,908175;
              30770100,911819; 30771000,910355; 30771900,909268; 30772800,
              909559; 30773700,908630; 30774600,910599; 30775500,898844;
              30776400,914461; 30777300,904651; 30778200,915404; 30779100,
              911512; 30780000,907986; 30780900,911760; 30781800,907949;
              30782700,908414; 30783600,904291; 30784500,915837; 30785400,
              908709; 30786300,908369; 30787200,910863; 30788100,910619;
              30789000,913489; 30789900,909668; 30790800,905649; 30791700,
              911220; 30792600,909513; 30793500,909696; 30794400,910039;
              30795300,912972; 30796200,909139; 30797100,909373; 30798000,
              910539; 30798900,911689; 30799800,908257; 30800700,909211;
              30801600,903026; 30802500,915334; 30803400,910699; 30804300,
              910145; 30805200,913410; 30806100,918210; 30807000,908803;
              30807900,901518; 30808800,906895; 30809700,913507; 30810600,
              905954; 30811500,912620; 30812400,904415; 30813300,912331;
              30814200,909121; 30815100,909953; 30816000,911005; 30816900,
              910311; 30817800,910597; 30818700,907942; 30819600,909653;
              30820500,909031; 30821400,907929; 30822300,910074; 30823200,
              914043; 30824100,905509; 30825000,911203; 30825900,908355;
              30826800,915409; 30827700,904642; 30828600,911345; 30829500,
              911000; 30830400,910084; 30831300,913044; 30832200,908071;
              30833100,911218; 30834000,907567; 30834900,911196; 30835800,
              905486; 30836700,913883; 30837600,912879; 30838500,908665;
              30839400,910092; 30840300,911341; 30841200,908811; 30842100,
              913347; 30843000,907413; 30843900,912818; 30844800,912178;
              30845700,906564; 30846600,909858; 30847500,912222; 30848400,
              910800; 30849300,911471; 30850200,910857; 30851100,909461;
              30852000,910389; 30852900,908306; 30853800,909538; 30854700,
              912366; 30855600,908670; 30856500,908529; 30857400,912628;
              30858300,908890; 30859200,907066; 30860100,912626; 30861000,
              908136; 30861900,908593; 30862800,910304; 30863700,911854;
              30864600,905906; 30865500,910895; 30866400,908357; 30867300,
              909955; 30868200,910114; 30869100,908739; 30870000,910457;
              30870900,911487; 30871800,911479; 30872700,911968; 30873600,
              909065; 30874500,907970; 30875400,913079; 30876300,907452;
              30877200,908824; 30878100,910579; 30879000,911219; 30879900,
              912349; 30880800,908627; 30881700,910734; 30882600,908012;
              30883500,910652; 30884400,909840; 30885300,912540; 30886200,
              906025; 30887100,912357; 30888000,911437; 30888900,908869;
              30889800,911200; 30890700,911981; 30891600,909028; 30892500,
              909355; 30893400,913093; 30894300,908080; 30895200,909333;
              30896100,911832; 30897000,922672; 30897900,901069; 30898800,
              913632; 30899700,907387; 30900600,905962; 30901500,913886;
              30902400,906322; 30903300,913673; 30904200,909374; 30905100,
              909154; 30906000,911078; 30906900,911226; 30907800,909244;
              30908700,907504; 30909600,911411; 30910500,908609; 30911400,
              909666; 30912300,904317; 30913200,916972; 30914100,907447;
              30915000,909474; 30915900,913598; 30916800,909123; 30917700,
              909262; 30918600,911336; 30919500,910681; 30920400,908585;
              30921300,910017; 30922200,911094; 30923100,908110; 30924000,
              910344; 30924900,908068; 30925800,906481; 30926700,910760;
              30927600,910835; 30928500,910665; 30929400,911538; 30930300,
              909263; 30931200,904900; 30932100,916856; 30933000,910661;
              30933900,911666; 30934800,910880; 30935700,909316; 30936600,
              909145; 30937500,913516; 30938400,904594; 30939300,912204;
              30940200,911578; 30941100,908896; 30942000,912178; 30942900,
              911628; 30943800,910431; 30944700,914030; 30945600,913675;
              30946500,907968; 30947400,911795; 30948300,909254; 30949200,
              907681; 30950100,897179; 30951000,915515; 30951900,906454;
              30952800,911078; 30953700,905135; 30954600,907936; 30955500,
              912760; 30956400,912841; 30957300,909793; 30958200,912026;
              30959100,912125; 30960000,905874; 30960900,908861; 30961800,
              917610; 30962700,916987; 30963600,906968; 30964500,911452;
              30965400,908816; 30966300,911088; 30967200,912809; 30968100,
              908209; 30969000,910788; 30969900,911967; 30970800,908015;
              30971700,913753; 30972600,909866; 30973500,911341; 30974400,
              906583; 30975300,911935; 30976200,911620; 30977100,908423;
              30978000,912067; 30978900,908826; 30979800,911482; 30980700,
              910138; 30981600,910751; 30982500,908783; 30983400,910915;
              30984300,907488; 30985200,913427; 30986100,909356; 30987000,
              909699; 30987900,908811; 30988800,911567; 30989700,911940;
              30990600,907108; 30991500,912249; 30992400,907876; 30993300,
              910709; 30994200,909982; 30995100,911118; 30996000,909634;
              30996900,909096; 30997800,910204; 30998700,908210; 30999600,
              911699; 31000500,910792; 31001400,907897; 31002300,912149;
              31003200,910144; 31004100,909852; 31005000,911338; 31005900,
              909787; 31006800,911237; 31007700,908994; 31008600,908414;
              31009500,908636; 31010400,916385; 31011300,905352; 31012200,
              914205; 31013100,906681; 31014000,910251; 31014900,910679;
              31015800,911190; 31016700,906377],
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{176,-66},{196,-46}})));
        Modelica.Blocks.Sources.RealExpression P_Ste_S(y=boi.port_b2.p)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{178,-92},{198,-72}})));
        Modelica.Blocks.Math.Gain             iQfue1(k=0.89)
          annotation (Placement(transformation(extent={{-140,0},{-120,20}})));
      equation
        connect(sou.ports[1],fwPum. port_a)
          annotation (Line(points={{-60,-20},{-40,-20}}, color={0,127,255},
            thickness=0.5));
        connect(fwPum.port_b, boi.port_a2) annotation (Line(points={{-20,-20},{0,
                -20},{0,16},{8,16}},
                                color={0,127,255},
            thickness=0.5));
        connect(TAmb.port, boi.heatPort) annotation (Line(points={{6,94},{18,94},
                {18,30}}, color={191,0,0}));
        connect(QFue.y,iQfue. u)
          annotation (Line(points={{-169,-116},{-154,-116}},
                                                           color={0,0,127}));
        connect(QWat.y,iQWat. u)
          annotation (Line(points={{-169,-166},{-154,-166}}, color={0,0,127}));
        connect(QFlue.y,iQFlue. u)
          annotation (Line(points={{-61,-118},{-42,-118}},
                                                         color={0,0,127}));
        connect(QLoss.y,iQLoss. u)
          annotation (Line(points={{-61,-168},{-42,-168}}, color={0,0,127}));
        connect(iQComb.u,Qcomb. y)
          annotation (Line(points={{92,-118},{67,-118}},color={0,0,127}));
        connect(iQevap.u,Qevap. y)
          annotation (Line(points={{92,-168},{67,-168}},  color={0,0,127}));
        connect(iQWat.y,perQWat. u1)
          annotation (Line(points={{-131,-166},{-118,-166}},color={0,0,127}));
        connect(iQfue.y,perQWat. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-154},{-118,-154}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQLoss.y,perQLoss. u1)
          annotation (Line(points={{-19,-168},{-8,-168}},
                                                        color={0,0,127}));
        connect(iQfue.y,perQLoss. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-156},{-8,-156}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQFlue.y,perQFlue. u1)
          annotation (Line(points={{-19,-118},{-6,-118}},
                                                      color={0,0,127}));
        connect(iQfue.y,perQFlue. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-130},{-6,-130}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQComb.y,perQcomb. u1)
          annotation (Line(points={{115,-118},{128,-118}},
                                                         color={0,0,127}));
        connect(iQfue.y,perQcomb. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-130},{128,-130}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQevap.y,perQeva. u1)
          annotation (Line(points={{115,-168},{128,-168}}, color={0,0,127}));
        connect(iQfue.y,perQeva. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-156},{128,-156}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(perQcomb.y,diffQcomb. u2) annotation (Line(points={{151,-124},{158,-124},
                {158,-210},{72,-210}}, color={0,0,127}));
        connect(perQeva.y,diffQevap. u2) annotation (Line(points={{151,-162},{150,-162},
                {150,-202},{14,-202},{14,-210}}, color={0,0,127}));
        connect(perQWat.y,diffQevap. u1) annotation (Line(points={{-95,-160},{-88,-160},
                {-88,-236},{26,-236},{26,-222},{14,-222}}, color={244,125,35}));
        connect(perQFlue.y,diffQcomb. u1) annotation (Line(points={{17,-124},{40,-124},
                {40,-236},{82,-236},{82,-222},{72,-222}},  color={244,125,35}));
        connect(boi.port_b1, senTem.port_b) annotation (Line(
            points={{28,24},{60,24}},
            color={0,0,0},
            thickness=0.5));
        connect(pro.ports[1], senTem1.port_b) annotation (Line(
            points={{-60,40},{-50,40}},
            color={0,0,0},
            thickness=0.5));
        connect(senTem1.port_a, boi.port_a1) annotation (Line(
            points={{-30,40},{-2,40},{-2,24},{8,24}},
            color={0,0,0},
            thickness=0.5));
        connect(boi.port_b2, senTem2.port_b) annotation (Line(
            points={{28,16},{40,16},{40,-20},{60,-20}},
            color={238,46,47},
            thickness=0.5));
        connect(pLR.y[1], boi.y) annotation (Line(
            points={{-99,80},{-20,80},{-20,52},{0,52},{0,28},{6,28}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(senTem.port_a, res.port_a)
          annotation (Line(points={{80,24},{100,24}}, color={0,127,255}));
        connect(res.port_b, exh.ports[1])
          annotation (Line(points={{120,24},{130,24}}, color={0,127,255}));
        connect(senTem2.port_a, res1.port_a)
          annotation (Line(points={{80,-20},{100,-20}}, color={0,127,255}));
        connect(res1.port_b, sin.ports[1])
          annotation (Line(points={{120,-20},{130,-20}}, color={0,127,255}));
        connect(iQfue1.y, fwPum.m_flow_in)
          annotation (Line(points={{-119,10},{-30,10},{-30,-8}}, color={0,0,127}));
        connect(mFloFw.y[1], iQfue1.u)
          annotation (Line(points={{-159,10},{-142,10}}, color={0,0,127}));
        connect(TAirIn.y[1], pro.T_in)
          annotation (Line(points={{-159,44},{-82,44}}, color={0,0,127}));
        connect(TFw.y[1], sou.T_in) annotation (Line(points={{-157,-32},{-90,-32},{-90,
                -16},{-82,-16}}, color={0,0,127}));
        annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Calibration/BoilerFurnace.mos"
              "Simulate and plot"),
          experiment(
            StartTime=30499200,
            StopTime=31016700,
            Tolerance=1e-06,
            __Dymola_Algorithm="Dassl"),            Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
                  {160,100}})),                                        Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}}),
              graphics={
                Rectangle(
                extent={{-200,-100},{160,-240}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid), Text(
                extent={{-198,-82},{-92,-102}},
                textColor={28,108,200},
                textString="Steady state verification (Heat balance)"),
                Rectangle(
                extent={{160,100},{300,-240}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid)}));
      end HeatBalanceBoilerTest;

      model InputOutputTest
        parameter String filnam = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTrain.mos");

             parameter String Inputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTrain.mos");
           parameter String Outputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTrainOut.mos");

        Modelica.Blocks.Sources.CombiTimeTable pLR(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(filnam),
          verboseRead=true,
          columns=2:5,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-42,16},{-22,36}})));

        Modelica.Blocks.Sources.CombiTimeTable pLR1(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={4},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-44,-18},{-24,2}})));
        Modelica.Blocks.Sources.CombiTimeTable pLR2(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Outputs),
          verboseRead=true,
          columns={4},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-42,-46},{-22,-26}})));
        annotation (                                 experiment(
            StartTime=8683200,
            StopTime=9114300,
            __Dymola_Algorithm="Dassl"));
      end InputOutputTest;

      model HeatBalanceBoilerForPaper
          extends Modelica.Icons.Example;

        // Medium declarations
        package MediumWat =
            Buildings.Media.Specialized.Water.TemperatureDependentDensity
          "Water medium - port_a (inlet)";
        package MediumSte = Buildings.Media.Steam (
          p_default=900000,
          T_default=273.15+180,
          h_default=2777170)
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

             parameter String Inputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTrain.mos");
           parameter String Outputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/BoilerTrainOut.mos");

        Fluid.Sources.Boundary_pT           sou(
          redeclare package Medium = MediumWat,
          p=100000,
          use_T_in=true,
          T=303.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
        Fluid.Movers.FlowControlled_m_flow fwPum(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m_flow_nominal,
          addPowerToMedium=false,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal=dp_nominal) "Feed water pump"
          annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
        Fluid.Sources.Boundary_pT com(
          redeclare package Medium = MediumAir,
          use_T_in=true,
          T=573.15,
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{-80,30},{-60,50}})));
        HeatTransfer.Sources.FixedTemperature           TAmb(T(displayUnit="K") = 305)
          "Ambient temperature in boiler room"
          annotation (Placement(transformation(extent={{80,60},{60,80}})));
        Fluid.Sources.Boundary_pT exh(
          redeclare package Medium = MediumAir,
          T(displayUnit="K"),
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{150,14},{130,34}})));

        Fluid.Sources.Boundary_pT           sin(
          redeclare package Medium = MediumSte,
          p(displayUnit="bar") = 895000,
          T=453.15,
          nPorts=1)
          "Sink"
          annotation (Placement(transformation(extent={{150,-30},{130,-10}})));
        BoilerPolynomialFurnaceHeatBalance boi(
          m1_flow_nominal=10,
          m2_flow_nominal=m_flow_nominal,
          show_T=true,
          energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
          p_start=900000,
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
          annotation (Placement(transformation(extent={{8,10},{28,30}})));

        Modelica.Blocks.Sources.RealExpression QFue(y=boi.QFue_flow)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{-190,-126},{-170,-106}})));
        Modelica.Blocks.Sources.RealExpression QWat(y=boi.QWat_flow)
          "Water heat flow rate"
          annotation (Placement(transformation(extent={{-190,-176},{-170,-156}})));
        Modelica.Blocks.Sources.RealExpression QLoss(y=boi.heatPort.Q_flow)
          "Loss from boiler casing"
          annotation (Placement(transformation(extent={{-82,-178},{-62,-158}})));
        Modelica.Blocks.Sources.RealExpression QFlue(y=boi.qExhFlo.y)
          "Heat losses into flue gas"
          annotation (Placement(transformation(extent={{-82,-128},{-62,-108}})));
        Modelica.Blocks.Sources.RealExpression Qcomb(y=((boi.port_b1.h_outflow - boi.port_a1.h_outflow)
              *boi.port_a1.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{46,-128},{66,-108}})));
        Modelica.Blocks.Sources.RealExpression Qevap(y=((boi.port_b2.h_outflow - boi.port_a2.h_outflow)
              *boi.port_b2.m_flow))
          "Q calcuated based on the mflow and enthalpy of the combustion side ports"
          annotation (Placement(transformation(extent={{46,-178},{66,-158}})));
        Modelica.Blocks.Continuous.Integrator iQfue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-152,-126},{-132,-106}})));
        Modelica.Blocks.Continuous.Integrator iQWat(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-152,-176},{-132,-156}})));
        Modelica.Blocks.Continuous.Integrator iQFlue(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-40,-128},{-20,-108}})));
        Modelica.Blocks.Continuous.Integrator iQLoss(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{-40,-178},{-20,-158}})));
        Modelica.Blocks.Continuous.Integrator iQComb(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{94,-128},{114,-108}})));
        Modelica.Blocks.Continuous.Integrator iQevap(k=1/(3600*1000), y_start=1)
          annotation (Placement(transformation(extent={{94,-178},{114,-158}})));
        Modelica.Blocks.Math.Division perQWat
          annotation (Placement(transformation(extent={{-116,-150},{-96,-170}})));
        Modelica.Blocks.Math.Division perQLoss
          annotation (Placement(transformation(extent={{-6,-152},{14,-172}})));
        Modelica.Blocks.Math.Division perQFlue
          annotation (Placement(transformation(extent={{-4,-134},{16,-114}})));
        Modelica.Blocks.Math.Division perQcomb
          annotation (Placement(transformation(extent={{130,-134},{150,-114}})));
        Modelica.Blocks.Math.Division perQeva
          annotation (Placement(transformation(extent={{130,-152},{150,-172}})));
        Modelica.Blocks.Math.Add diffQevap(k1=-1)
          annotation (Placement(transformation(extent={{12,-206},{-8,-226}})));
        Modelica.Blocks.Math.Add diffQcomb(k1=-1)
          annotation (Placement(transformation(extent={{70,-206},{50,-226}})));
        Fluid.Sensors.TemperatureTwoPort TExh(
          redeclare package Medium = MediumAir,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K")) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={70,24})));
        Fluid.Sensors.TemperatureTwoPort
                                       senTem1(
          redeclare package Medium = MediumAir,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K"))
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=180,
              origin={-40,40})));
        Fluid.Sensors.TemperatureTwoPort TSte(
          redeclare package Medium = MediumSte,
          m_flow_nominal=m_flow_nominal,
          tau=30,
          T_start(displayUnit="K")) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={70,-20})));
        Modelica.Blocks.Sources.RealExpression Qboi_S(y=QWat.y - QLoss.y)
          "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{176,28},{196,48}})));
        Modelica.Blocks.Sources.RealExpression mFueFlo_S(y=boi.mFue_flow)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,8},{196,28}})));
        Modelica.Blocks.Sources.RealExpression T_exh_S(y=TExh.T)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,-12},{196,8}})));
        Modelica.Blocks.Sources.RealExpression T_Ste_S(y=TSte.T)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{176,-32},{196,-12}})));
        Modelica.Blocks.Sources.CombiTimeTable plr(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={5},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-120,68},{-100,88}})));
        Modelica.Blocks.Sources.CombiTimeTable TAirIn(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={2},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-120,34},{-100,54}})));
        Modelica.Blocks.Sources.CombiTimeTable mFloFw(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={3},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-120,0},{-100,20}})));
        Modelica.Blocks.Sources.CombiTimeTable Tfw(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={4},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-120,-26},{-100,-6}})));
        Modelica.Blocks.Sources.CombiTimeTable MeaData(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Outputs),
          verboseRead=true,
          columns=2:5,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{176,62},{196,82}})));
        Fluid.FixedResistances.PressureDrop res(
          redeclare package Medium = MediumAir,
          m_flow_nominal=10,
          dp_nominal=50000)
          annotation (Placement(transformation(extent={{100,14},{120,34}})));
        Fluid.FixedResistances.PressureDrop res1(
          redeclare package Medium = MediumSte,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=10000)
          annotation (Placement(transformation(extent={{100,-30},{120,-10}})));
        Modelica.Blocks.Sources.CombiTimeTable stePre(
          tableOnFile=false,
          table=[8683200,892316; 8684100,891685; 8685000,890902; 8685900,891196; 8686800,
              895785; 8687700,898587; 8688600,900129; 8689500,896132; 8690400,891488;
              8691300,885147; 8692200,889115; 8693100,896804; 8694000,893937; 8694900,
              894953; 8695800,892835; 8696700,891911; 8697600,893801; 8698500,897013;
              8699400,890919; 8700300,894550; 8701200,894498; 8702100,891167; 8703000,
              887755; 8703900,891254; 8704800,891116; 8705700,895712; 8706600,891681;
              8707500,890279; 8708400,893290; 8709300,895085; 8710200,893762; 8711100,
              893850; 8712000,891206; 8712900,889010; 8713800,891960; 8714700,893310;
              8715600,894229; 8716500,894933; 8717400,892706; 8718300,890688; 8719200,
              892198; 8720100,888249; 8721000,890422; 8721900,895005; 8722800,894123;
              8723700,893963; 8724600,892238; 8725500,891018; 8726400,893750; 8727300,
              892906; 8728200,891890; 8729100,892280; 8730000,890381; 8730900,889701;
              8731800,892948; 8732700,890220; 8733600,892233; 8734500,891652; 8735400,
              888630; 8736300,890460; 8737200,890518; 8738100,890228; 8739000,895814;
              8739900,900284; 8740800,900628; 8741700,898741; 8742600,898387; 8743500,
              896780; 8744400,895804; 8745300,896948; 8746200,896786; 8747100,897274;
              8748000,892660; 8748900,897156; 8749800,896337; 8750700,897688; 8751600,
              892595; 8752500,897295; 8753400,898425; 8754300,894919; 8755200,899118;
              8756100,896767; 8757000,897065; 8757900,894714; 8758800,889521; 8759700,
              901560; 8760600,900437; 8761500,901380; 8762400,898418; 8763300,899716;
              8764200,897777; 8765100,900712; 8766000,898616; 8766900,898191; 8767800,
              897273; 8768700,900192; 8769600,902553; 8770500,901377; 8771400,901969;
              8772300,899335; 8773200,900886; 8774100,899350; 8775000,900286; 8775900,
              898729; 8776800,902840; 8777700,899746; 8778600,900780; 8779500,901832;
              8780400,901287; 8781300,901840; 8782200,901655; 8783100,900987; 8784000,
              899158; 8784900,899444; 8785800,897361; 8786700,897411; 8787600,896690;
              8788500,894345; 8789400,897702; 8790300,899003; 8791200,899007; 8792100,
              897467; 8793000,898909; 8793900,896982; 8794800,899785; 8795700,898488;
              8796600,897608; 8797500,901257; 8798400,898667; 8799300,897137; 8800200,
              899305; 8801100,899503; 8802000,899118; 8802900,900084; 8803800,898564;
              8804700,899755; 8805600,902665; 8806500,898058; 8807400,901453; 8808300,
              900676; 8809200,897994; 8810100,902194; 8811000,898462; 8811900,897836;
              8812800,902642; 8813700,902053; 8814600,898060; 8815500,898509; 8816400,
              901198; 8817300,902410; 8818200,900042; 8819100,901264; 8820000,901855;
              8820900,903825; 8821800,898006; 8822700,902393; 8823600,897728; 8824500,
              897881; 8825400,899065; 8826300,899409; 8827200,900042; 8828100,899106;
              8829000,894843; 8829900,897442; 8830800,891668; 8831700,899983; 8832600,
              900388; 8833500,898621; 8834400,894455; 8835300,899183; 8836200,902382;
              8837100,901714; 8838000,892222; 8838900,897753; 8839800,895468; 8840700,
              896568; 8841600,900810; 8842500,898862; 8843400,902657; 8844300,907208;
              8845200,901789; 8846100,905161; 8847000,904069; 8847900,919049; 8848800,
              933409; 8849700,933443; 8850600,937552; 8851500,917675; 8852400,914861;
              8853300,927014; 8854200,924750; 8855100,920426; 8856000,909951; 8856900,
              913132; 8857800,920944; 8858700,905232; 8859600,904585; 8860500,903512;
              8861400,903256; 8862300,899210; 8863200,899369; 8864100,897609; 8865000,
              899453; 8865900,899356; 8866800,898987; 8867700,898571; 8868600,897670;
              8869500,898591; 8870400,898527; 8871300,897947; 8872200,892801; 8873100,
              891317; 8874000,898983; 8874900,892916; 8875800,897873; 8876700,896933;
              8877600,895837; 8878500,896179; 8879400,900107; 8880300,897665; 8881200,
              898789; 8882100,906752; 8883000,898733; 8883900,896273; 8884800,903821;
              8885700,896866; 8886600,898431; 8887500,901278; 8888400,882846; 8889300,
              895956; 8890200,902357; 8891100,901307; 8892000,899826; 8892900,900040;
              8893800,896459; 8894700,897382; 8895600,899511; 8896500,896183; 8897400,
              903493; 8898300,899879; 8899200,899473; 8900100,902817; 8901000,897306;
              8901900,901992; 8902800,900536; 8903700,903393; 8904600,901198; 8905500,
              902132; 8906400,899723; 8907300,897768; 8908200,898384; 8909100,894115;
              8910000,894856; 8910900,892514; 8911800,901170; 8912700,896023; 8913600,
              901879; 8914500,897391; 8915400,899412; 8916300,904120; 8917200,899862;
              8918100,894087; 8919000,895448; 8919900,895685; 8920800,909290; 8921700,
              905152; 8922600,895080; 8923500,900401; 8924400,902049; 8925300,897232;
              8926200,895718; 8927100,900294; 8928000,901106; 8928900,901236; 8929800,
              902719; 8930700,900085; 8931600,902612; 8932500,901231; 8933400,902378;
              8934300,884259; 8935200,888231; 8936100,899815; 8937000,893617; 8937900,
              896661; 8938800,902966; 8939700,908549; 8940600,901597; 8941500,894364;
              8942400,895806; 8943300,901993; 8944200,893605; 8945100,887119; 8946000,
              894041; 8946900,898020; 8947800,898204; 8948700,901734; 8949600,906864;
              8950500,896506; 8951400,895502; 8952300,890127; 8953200,902039; 8954100,
              902657; 8955000,881041; 8955900,881701; 8956800,907420; 8957700,901602;
              8958600,900108; 8959500,895231; 8960400,896903; 8961300,899413; 8962200,
              895501; 8963100,895750; 8964000,898890; 8964900,897228; 8965800,896884;
              8966700,896233; 8967600,897827; 8968500,896559; 8969400,896825; 8970300,
              896490; 8971200,898713; 8972100,896536; 8973000,895528; 8973900,896828;
              8974800,897874; 8975700,897805; 8976600,900293; 8977500,895008; 8978400,
              894395; 8979300,897377; 8980200,896797; 8981100,901117; 8982000,901607;
              8982900,896163; 8983800,890473; 8984700,895937; 8985600,900544; 8986500,
              899973; 8987400,901876; 8988300,899390; 8989200,896199; 8990100,898284;
              8991000,899128; 8991900,901099; 8992800,898263; 8993700,896992; 8994600,
              900369; 8995500,899205; 8996400,901288; 8997300,903209; 8998200,899505;
              8999100,898207; 9000000,889801; 9000900,893753; 9001800,900226; 9002700,
              899721; 9003600,902705; 9004500,904361; 9005400,902308; 9006300,889626;
              9007200,885060; 9008100,894836; 9009000,903038; 9009900,894932; 9010800,
              893072; 9011700,906769; 9012600,904800; 9013500,901890; 9014400,900425;
              9015300,907985; 9016200,911117; 9017100,913138; 9018000,921595; 9018900,
              938608; 9019800,952207; 9020700,923140; 9021600,900445; 9022500,905692;
              9023400,908128; 9024300,889486; 9025200,897555; 9026100,896083; 9027000,
              888324; 9027900,889318; 9028800,889185; 9029700,898645; 9030600,899778;
              9031500,897199; 9032400,898766; 9033300,898734; 9034200,896020; 9035100,
              891446; 9036000,899859; 9036900,899089; 9037800,897368; 9038700,898364;
              9039600,896495; 9040500,894885; 9041400,894342; 9042300,899818; 9043200,
              899519; 9044100,896962; 9045000,895953; 9045900,894465; 9046800,895224;
              9047700,894825; 9048600,895060; 9049500,897563; 9050400,894841; 9051300,
              894129; 9052200,897354; 9053100,897693; 9054000,899312; 9054900,895058;
              9055800,894160; 9056700,896299; 9057600,896457; 9058500,897513; 9059400,
              898402; 9060300,897171; 9061200,895457; 9062100,897859; 9063000,897867;
              9063900,897006; 9064800,899385; 9065700,895948; 9066600,894770; 9067500,
              897424; 9068400,894896; 9069300,892490; 9070200,890507; 9071100,894937;
              9072000,897981; 9072900,898133; 9073800,901646; 9074700,901137; 9075600,
              896753; 9076500,890813; 9077400,896331; 9078300,896888; 9079200,904427;
              9080100,901888; 9081000,893080; 9081900,891112; 9082800,895399; 9083700,
              897208; 9084600,900491; 9085500,894977; 9086400,893689; 9087300,893570;
              9088200,890527; 9089100,894471; 9090000,891748; 9090900,894354; 9091800,
              896214; 9092700,894799; 9093600,890308; 9094500,893537; 9095400,893195;
              9096300,890145; 9097200,891324; 9098100,886895; 9099000,888348; 9099900,
              902638; 9100800,896894; 9101700,894346; 9102600,893451; 9103500,893956;
              9104400,904706; 9105300,900887; 9106200,895223; 9107100,892098; 9108000,
              898905; 9108900,896420; 9109800,896890; 9110700,900077; 9111600,901466;
              9112500,898862; 9113400,900364; 9114300,899451],
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{176,-66},{196,-46}})));
        Modelica.Blocks.Sources.RealExpression P_Ste_S(y=boi.port_b2.p)
          "Fuel mass flow rate"
          annotation (Placement(transformation(extent={{178,-92},{198,-72}})));
      equation
        connect(sou.ports[1],fwPum. port_a)
          annotation (Line(points={{-60,-20},{-40,-20}}, color={0,127,255},
            thickness=0.5));
        connect(fwPum.port_b, boi.port_a2) annotation (Line(points={{-20,-20},{0,
                -20},{0,16},{8,16}},
                                color={0,127,255},
            thickness=0.5));
        connect(TAmb.port, boi.heatPort) annotation (Line(points={{60,70},{18,
                70},{18,30}},
                          color={191,0,0}));
        connect(QFue.y,iQfue. u)
          annotation (Line(points={{-169,-116},{-154,-116}},
                                                           color={0,0,127}));
        connect(QWat.y,iQWat. u)
          annotation (Line(points={{-169,-166},{-154,-166}}, color={0,0,127}));
        connect(QFlue.y,iQFlue. u)
          annotation (Line(points={{-61,-118},{-42,-118}},
                                                         color={0,0,127}));
        connect(QLoss.y,iQLoss. u)
          annotation (Line(points={{-61,-168},{-42,-168}}, color={0,0,127}));
        connect(iQComb.u,Qcomb. y)
          annotation (Line(points={{92,-118},{67,-118}},color={0,0,127}));
        connect(iQevap.u,Qevap. y)
          annotation (Line(points={{92,-168},{67,-168}},  color={0,0,127}));
        connect(iQWat.y,perQWat. u1)
          annotation (Line(points={{-131,-166},{-118,-166}},color={0,0,127}));
        connect(iQfue.y,perQWat. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-154},{-118,-154}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQLoss.y,perQLoss. u1)
          annotation (Line(points={{-19,-168},{-8,-168}},
                                                        color={0,0,127}));
        connect(iQfue.y,perQLoss. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-156},{-8,-156}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQFlue.y,perQFlue. u1)
          annotation (Line(points={{-19,-118},{-6,-118}},
                                                      color={0,0,127}));
        connect(iQfue.y,perQFlue. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{-14,-142},{-14,-130},{-6,-130}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQComb.y,perQcomb. u1)
          annotation (Line(points={{115,-118},{128,-118}},
                                                         color={0,0,127}));
        connect(iQfue.y,perQcomb. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-130},{128,-130}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(iQevap.y,perQeva. u1)
          annotation (Line(points={{115,-168},{128,-168}}, color={0,0,127}));
        connect(iQfue.y,perQeva. u2) annotation (Line(
            points={{-131,-116},{-124,-116},{-124,-142},{120,-142},{120,-156},{128,-156}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(perQcomb.y,diffQcomb. u2) annotation (Line(points={{151,-124},{158,-124},
                {158,-210},{72,-210}}, color={0,0,127}));
        connect(perQeva.y,diffQevap. u2) annotation (Line(points={{151,-162},{150,-162},
                {150,-202},{14,-202},{14,-210}}, color={0,0,127}));
        connect(perQWat.y,diffQevap. u1) annotation (Line(points={{-95,-160},{-88,-160},
                {-88,-236},{26,-236},{26,-222},{14,-222}}, color={244,125,35}));
        connect(perQFlue.y,diffQcomb. u1) annotation (Line(points={{17,-124},{40,-124},
                {40,-236},{82,-236},{82,-222},{72,-222}},  color={244,125,35}));
        connect(boi.port_b1, TExh.port_b) annotation (Line(
            points={{28,24},{60,24}},
            color={0,0,0},
            thickness=0.5));
        connect(com.ports[1], senTem1.port_b) annotation (Line(
            points={{-60,40},{-50,40}},
            color={0,0,0},
            thickness=0.5));
        connect(senTem1.port_a, boi.port_a1) annotation (Line(
            points={{-30,40},{-2,40},{-2,24},{8,24}},
            color={0,0,0},
            thickness=0.5));
        connect(boi.port_b2, TSte.port_b) annotation (Line(
            points={{28,16},{40,16},{40,-20},{60,-20}},
            color={238,46,47},
            thickness=0.5));
        connect(plr.y[1], boi.y) annotation (Line(
            points={{-99,78},{-20,78},{-20,52},{0,52},{0,28},{6,28}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(TExh.port_a, res.port_a) annotation (Line(
            points={{80,24},{100,24}},
            color={0,0,0},
            thickness=0.5));
        connect(res.port_b, exh.ports[1])
          annotation (Line(points={{120,24},{130,24}}, color={0,0,0},
            thickness=0.5));
        connect(TSte.port_a, res1.port_a) annotation (Line(
            points={{80,-20},{100,-20}},
            color={238,46,47},
            thickness=0.5));
        connect(res1.port_b, sin.ports[1])
          annotation (Line(points={{120,-20},{130,-20}}, color={238,46,47},
            thickness=0.5));
        connect(TAirIn.y[1], com.T_in) annotation (Line(
            points={{-99,44},{-82,44}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(Tfw.y[1], sou.T_in) annotation (Line(
            points={{-99,-16},{-82,-16}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(mFloFw.y[1], fwPum.m_flow_in) annotation (Line(
            points={{-99,10},{-30,10},{-30,-8}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Calibration/BoilerFurnace.mos"
              "Simulate and plot"),
          experiment(
            StartTime=8683200,
            StopTime=9114336,
            Tolerance=1e-06,
            __Dymola_Algorithm="Dassl"),            Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
                  {160,100}})),                                        Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}}),
              graphics={
                Rectangle(
                extent={{-200,-100},{160,-240}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid), Text(
                extent={{-198,-82},{-92,-102}},
                textColor={28,108,200},
                textString="Steady state verification (Heat balance)"),
                Rectangle(
                extent={{160,100},{300,-240}},
                lineColor={28,108,200},
                fillColor={174,179,179},
                fillPattern=FillPattern.Solid)}));
      end HeatBalanceBoilerForPaper;

      package Optimization

        function OptimizationLastRunModel
        /* Automatically generated at Mon Apr  3 20:33:38 2023 */

        /* This function needs Optimization library 2.2.3 or higher */

            input Boolean interactive=true annotation(Dialog(tab="Advanced"));
            output Boolean runOK annotation(Dialog(tab="Advanced", group="Output"));
        protected
            Optimization.Internal.Version.V22.ModelOptimizationSetup setup=
                Optimization.Internal.Version.V22.ModelOptimizationSetup(modelName="Buildings.GEDCalibration.CUBoulder.Components.Validation.Calibration.EconomizerTrain",
            plotScript="",
            saveSetup=true,
            saveSetupFilename="OptimizationLastRunModel.mo",
            convertSetup=false,
            askForTunerReUse=true,
            tuner=Optimization.Internal.Version.Current.Tuner(
                UseTunerMatrixForDiscreteValues=false,
                tunerParameters={Optimization.Internal.Version.Current.TunerParameter(
                  name="Q_flow_nominal",
                  active=true,
                  Value=1324000,
                  scaleToBounds=false,
                  min=1100000,
                  max=2400000,
                  equidistant=0,
                  discreteValues=fill(0.0, 0),
                  unit="W"),Optimization.Internal.Version.Current.TunerParameter(
                  name="r_nominal",
                  active=true,
                  Value=0.666667,
                  min=0,
                  max=1,
                  discreteValues=fill(0, 0),
                  unit=""),Optimization.Internal.Version.Current.TunerParameter(
                  name="T_a2_nominal",
                  active=true,
                  Value=250,
                  min=200,
                  max=300,
                  discreteValues=fill(0, 0),
                  unit="K"),Optimization.Internal.Version.Current.TunerParameter(
                  name="T_a1_nominal",
                  active=true,
                  Value=577,
                  min=450,
                  max=650,
                  discreteValues=fill(0, 0),
                  unit="K")},
                tunerMatrix=fill(
                  0.0,
                  0,
                  4)),
            criteria={Optimization.Internal.Version.Current.Criterion(
                name="diffFW.y",
                active=true,
                usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                demand=1.0,
                unit=""),Optimization.Internal.Version.Current.Criterion(
                name="diffDfT.y",
                active=true,
                usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                demand=1,
                unit="")},
            preferences=Optimization.Internal.Version.Current.Preferences(
                optimizationOptions=
                  Optimization.Internal.Version.Current.OptimizationOptions(
                  method=Optimization.Internal.Version.Current.Types.OptimizationMethod.simsa,
                  ObjectiveFunctionType=Optimization.Internal.Version.Current.Types.ObjectiveFunctionType.Max,
                  OptTol=0.00001,
                  maxEval=1000,
                  GridBlock=50,
                  evalBestFinal=false,
                  saveBest=true,
                  saveHistory=true,
                  listFilename="OptimizationLog.log",
                  listOn=true,
                  listOnline=true,
                  listIncrement=100,
                  numberOfShownDigits=3,
                  onPlace=true,
                  listTuners=true,
                  GaPopSize=10,
                  GaNGen=100),
                simulationOptions=
                  Optimization.Internal.Version.Current.SimulationOptions(
                  startTime=0.0,
                  stopTime=1.0,
                  outputInterval=0.0,
                  numberOfIntervals=500,
                  integrationMethod=Optimization.Internal.Version.Current.Types.IntegrationMethod.Dassl,
                  integrationTolerance=0.001,
                  fixedStepSize=0.0,
                  autoLoadResults=true,
                  useDsFinal=true,
                  translateModel=false,
                  setCriteriaSimulationFailed=true,
                  CriteriaSimulationFailedValue=1000000.0,
                  simulationMode=Optimization.Internal.Version.Current.Types.SimulationMode.Single,
                  parallelizationMode=Optimization.Internal.Version.Current.Types.ParallelizationMode.None,
                  numberOfThreads=0,
                  copyFiles=fill("", 0)),
                sensitivityOptions=
                  Optimization.Internal.Version.Current.SensitivityOptions(
                  TypeOfSensitivityComputation=Optimization.Internal.Version.Current.Types.SensitivityMethod.ExternalDifferencesSymmetric,
                  automaticSensitivityTolerance=true,
                  sensitivityTolerance=1E-06)));

        algorithm
            runOK := Optimization.Tasks.ModelOptimization.run22(setup, interactive);
            annotation(__Dymola_interactive=true, __Dymola_DymolaStoredErrors(
                thetext="function OptimizationLastRunModel
/* Automatically generated at Mon Apr  3 20:33:38 2023 */

/* This function needs Optimization library 2.2.3 or higher */

    input Boolean interactive=true annotation(Dialog(tab=\"Advanced\"));
    output Boolean runOK annotation(Dialog(tab=\"Advanced\", group=\"Output\"));
protected 
    Optimization.Internal.Version.V22.ModelOptimizationSetup setup=
        Optimization.Internal.Version.V22.ModelOptimizationSetup(

algorithm 
    runOK := Optimization.Tasks.ModelOptimization.run22(setup, interactive);
    annotation(__Dymola_interactive=true);
end OptimizationLastRunModel;
"));
        end OptimizationLastRunModel;

        record EcoOpt
          import Optimization;
          extends
            Optimization.Internal.Version.Current.ModelOptimizationSetup(
            modelName=
                "Buildings.GEDCalibration.CUBoulder.Components.Validation.Calibration.EconomizerTrain",
            plotScript="",
            saveSetup=true,
            saveSetupFilename="OptimizationLastRunModel.mo",
            convertSetup=false,
            askForTunerReUse=true,
            tuner=Optimization.Internal.Version.Current.Tuner(
                            UseTunerMatrixForDiscreteValues=false,
                            tunerParameters={
                  Optimization.Internal.Version.Current.TunerParameter(
                              name="Q_flow_nominal",
                              active=true,
                              Value=1324000,
                              scaleToBounds=false,
                              min=1100000,
                              max=2400000,
                              equidistant=0,
                              discreteValues=fill(0.0, 0),
                              unit="W"),
                  Optimization.Internal.Version.Current.TunerParameter(
                              name="r_nominal",
                              active=true,
                              Value=0.666667,
                              min=0,
                              max=1,
                              discreteValues=fill(0, 0),
                              unit=""),
                  Optimization.Internal.Version.Current.TunerParameter(
                              name="T_a2_nominal",
                              active=true,
                              Value=250,
                              min=200,
                              max=300,
                              discreteValues=fill(0, 0),
                              unit="K"),
                  Optimization.Internal.Version.Current.TunerParameter(
                              name="T_a1_nominal",
                              active=true,
                              Value=577,
                              min=450,
                              max=650,
                              discreteValues=fill(0, 0),
                              unit="K")},
                            tunerMatrix=fill(
                              0.0,
                              0,
                              4)),
            criteria={Optimization.Internal.Version.Current.Criterion(
                            name="diffFW.y",
                            active=true,
                            usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                            demand=1.0,
                            unit=""),
                Optimization.Internal.Version.Current.Criterion(
                            name="diffDfT.y",
                            active=true,
                            usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                            demand=1,
                            unit="")},
            preferences=Optimization.Internal.Version.Current.Preferences(
                            optimizationOptions=
                  Optimization.Internal.Version.Current.OptimizationOptions(
                              method=Optimization.Internal.Version.Current.Types.OptimizationMethod.simsa,
                              ObjectiveFunctionType=Optimization.Internal.Version.Current.Types.ObjectiveFunctionType.Max,
                              OptTol=0.00001,
                              maxEval=1000,
                              GridBlock=50,
                              evalBestFinal=false,
                              saveBest=true,
                              saveHistory=true,
                              listFilename="OptimizationLog.log",
                              listOn=true,
                              listOnline=true,
                              listIncrement=100,
                              numberOfShownDigits=3,
                              onPlace=true,
                              listTuners=true,
                              GaPopSize=10,
                              GaNGen=100),
                            simulationOptions=
                  Optimization.Internal.Version.Current.SimulationOptions(
                              startTime=0.0,
                              stopTime=1.0,
                              outputInterval=0.0,
                              numberOfIntervals=500,
                              integrationMethod=Optimization.Internal.Version.Current.Types.IntegrationMethod.Dassl,
                              integrationTolerance=0.001,
                              fixedStepSize=0.0,
                              autoLoadResults=true,
                              useDsFinal=true,
                              translateModel=false,
                              setCriteriaSimulationFailed=true,
                              CriteriaSimulationFailedValue=1000000.0,
                              simulationMode=Optimization.Internal.Version.Current.Types.SimulationMode.Single,
                              parallelizationMode=Optimization.Internal.Version.Current.Types.ParallelizationMode.None,
                              numberOfThreads=0,
                              copyFiles=fill("", 0)),
                            sensitivityOptions=
                  Optimization.Internal.Version.Current.SensitivityOptions(
                              TypeOfSensitivityComputation=Optimization.Internal.Version.Current.Types.SensitivityMethod.ExternalDifferencesSymmetric,
                              automaticSensitivityTolerance=true,
                              sensitivityTolerance=1E-06)));

        end EcoOpt;

        record EconomizerOptimization
          extends
            Optimization.Internal.Version.Current.ModelOptimizationSetup(
            modelName=
                "Buildings.GEDCalibration.CUBoulder.Components.Validation.Calibration.EconomizerTrain",
            plotScript="",
            saveSetup=true,
            saveSetupFilename="OptimizationLastRunModel.mo",
            convertSetup=false,
            askForTunerReUse=false,
            tuner=Optimization.Internal.Version.Current.Tuner(
                      UseTunerMatrixForDiscreteValues=false,
                      tunerParameters={
                  Optimization.Internal.Version.Current.TunerParameter(
                        name="Q_flow_nominal",
                        active=true,
                        Value=1324000,
                        scaleToBounds=false,
                        min=1024000,
                        max=2240000,
                        equidistant=0,
                        discreteValues=fill(0.0, 0),
                        unit="W")},
                      tunerMatrix=fill(
                        0.0,
                        0,
                        1)),
            criteria={Optimization.Internal.Version.Current.Criterion(
                      name="diffFW.y",
                      active=true,
                      usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                      demand=1.0,
                      unit=""),Optimization.Internal.Version.Current.Criterion(
                      name="diffDfT.y",
                      active=true,
                      usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                      demand=1,
                      unit="")},
            preferences=Optimization.Internal.Version.Current.Preferences(
                      optimizationOptions=
                  Optimization.Internal.Version.Current.OptimizationOptions(
                        method=Optimization.Internal.Version.Current.Types.OptimizationMethod.simsa,
                        ObjectiveFunctionType=Optimization.Internal.Version.Current.Types.ObjectiveFunctionType.Max,
                        OptTol=0.0001,
                        maxEval=1000,
                        GridBlock=50,
                        evalBestFinal=false,
                        saveBest=true,
                        saveHistory=true,
                        listFilename="OptimizationLog.log",
                        listOn=true,
                        listOnline=true,
                        listIncrement=100,
                        numberOfShownDigits=3,
                        onPlace=true,
                        listTuners=true,
                        GaPopSize=10,
                        GaNGen=100),
                      simulationOptions=
                  Optimization.Internal.Version.Current.SimulationOptions(
                        startTime=0.0,
                        stopTime=1.0,
                        outputInterval=0.0,
                        numberOfIntervals=500,
                        integrationMethod=Optimization.Internal.Version.Current.Types.IntegrationMethod.Dassl,
                        integrationTolerance=0.001,
                        fixedStepSize=0.0,
                        autoLoadResults=true,
                        useDsFinal=true,
                        translateModel=false,
                        setCriteriaSimulationFailed=true,
                        CriteriaSimulationFailedValue=1000000.0,
                        simulationMode=Optimization.Internal.Version.Current.Types.SimulationMode.Single,
                        parallelizationMode=Optimization.Internal.Version.Current.Types.ParallelizationMode.None,
                        numberOfThreads=0,
                        copyFiles=fill("", 0)),
                      sensitivityOptions=
                  Optimization.Internal.Version.Current.SensitivityOptions(
                        TypeOfSensitivityComputation=Optimization.Internal.Version.Current.Types.SensitivityMethod.ExternalDifferencesSymmetric,
                        automaticSensitivityTolerance=true,
                        sensitivityTolerance=1E-06)));

        end EconomizerOptimization;
      end Optimization;

      function OptimizationLastRunModel
        import Optimization;
      /* Automatically generated at Mon Apr  3 20:33:38 2023 */

      /* This function needs Optimization library 2.2.3 or higher */

          input Boolean interactive=true annotation(Dialog(tab="Advanced"));
          output Boolean runOK annotation(Dialog(tab="Advanced", group="Output"));
      protected
          Optimization.Internal.Version.V22.ModelOptimizationSetup setup=
            Optimization.Internal.Version.V22.ModelOptimizationSetup(
                  modelName=
              "Buildings.GEDCalibration.CUBoulder.Components.Validation.Calibration.EconomizerTrain",
                  plotScript="",
                  saveSetup=true,
                  saveSetupFilename="OptimizationLastRunModel.mo",
                  convertSetup=false,
                  askForTunerReUse=true,
                  tuner=Optimization.Internal.Version.Current.Tuner(
                    UseTunerMatrixForDiscreteValues=false,
                    tunerParameters={
                Optimization.Internal.Version.Current.TunerParameter(
                      name="Q_flow_nominal",
                      active=true,
                      Value=1324000,
                      scaleToBounds=false,
                      min=1100000,
                      max=2400000,
                      equidistant=0,
                      discreteValues=fill(0.0, 0),
                      unit="W"),
                Optimization.Internal.Version.Current.TunerParameter(
                      name="r_nominal",
                      active=true,
                      Value=0.666667,
                      min=0,
                      max=1,
                      discreteValues=fill(0, 0),
                      unit=""),
                Optimization.Internal.Version.Current.TunerParameter(
                      name="T_a2_nominal",
                      active=true,
                      Value=250,
                      min=200,
                      max=300,
                      discreteValues=fill(0, 0),
                      unit="K"),
                Optimization.Internal.Version.Current.TunerParameter(
                      name="T_a1_nominal",
                      active=true,
                      Value=577,
                      min=450,
                      max=650,
                      discreteValues=fill(0, 0),
                      unit="K")},
                    tunerMatrix=fill(
                      0.0,
                      0,
                      4)),
                  criteria={Optimization.Internal.Version.Current.Criterion(
                    name="diffFW.y",
                    active=true,
                    usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                    demand=1.0,
                    unit=""),Optimization.Internal.Version.Current.Criterion(
                    name="diffDfT.y",
                    active=true,
                    usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                    demand=1,
                    unit="")},
                  preferences=Optimization.Internal.Version.Current.Preferences(
                    optimizationOptions=
                Optimization.Internal.Version.Current.OptimizationOptions(
                      method=Optimization.Internal.Version.Current.Types.OptimizationMethod.simsa,
                      ObjectiveFunctionType=Optimization.Internal.Version.Current.Types.ObjectiveFunctionType.Max,
                      OptTol=0.00001,
                      maxEval=1000,
                      GridBlock=50,
                      evalBestFinal=false,
                      saveBest=true,
                      saveHistory=true,
                      listFilename="OptimizationLog.log",
                      listOn=true,
                      listOnline=true,
                      listIncrement=100,
                      numberOfShownDigits=3,
                      onPlace=true,
                      listTuners=true,
                      GaPopSize=10,
                      GaNGen=100),
                    simulationOptions=
                Optimization.Internal.Version.Current.SimulationOptions(
                      startTime=0.0,
                      stopTime=1.0,
                      outputInterval=0.0,
                      numberOfIntervals=500,
                      integrationMethod=Optimization.Internal.Version.Current.Types.IntegrationMethod.Dassl,
                      integrationTolerance=0.001,
                      fixedStepSize=0.0,
                      autoLoadResults=true,
                      useDsFinal=true,
                      translateModel=false,
                      setCriteriaSimulationFailed=true,
                      CriteriaSimulationFailedValue=1000000.0,
                      simulationMode=Optimization.Internal.Version.Current.Types.SimulationMode.Single,
                      parallelizationMode=Optimization.Internal.Version.Current.Types.ParallelizationMode.None,
                      numberOfThreads=0,
                      copyFiles=fill("", 0)),
                    sensitivityOptions=
                Optimization.Internal.Version.Current.SensitivityOptions(
                      TypeOfSensitivityComputation=Optimization.Internal.Version.Current.Types.SensitivityMethod.ExternalDifferencesSymmetric,
                      automaticSensitivityTolerance=true,
                      sensitivityTolerance=1E-06)));
      algorithm
          runOK := Optimization.Tasks.ModelOptimization.run22(setup, interactive);
          annotation(__Dymola_interactive=true);
      end OptimizationLastRunModel;
    end Calibration;

    package SteamCoilHeater

    end SteamCoilHeater;

    package ExhaustEconomizer
      model EconomizerTrain
          extends Modelica.Icons.Example;

           parameter String Inputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/EcoTrainIn.mos");
           parameter String Outputs = ("modelica://Buildings/Resources/Data/Experimental/GEDCalibration/EcoTrainOut.mos");

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
        parameter Modelica.Units.SI.Power Q_flow_nominal = 1324000 "Nominal power";
        parameter Modelica.Units.SI.SpecificEnthalpy dh_nominal=
          MediumSte.specificEnthalpy(
            MediumSte.setState_pTX(p=p_nominal, T=T_nominal, X=MediumSte.X_default))
          "Nominal change in enthalpy";
        parameter Modelica.Units.SI.MassFlowRate m1_flow_nominal=3.43
          "Nominal mass flow rate";

        parameter Modelica.Units.SI.MassFlowRate m2_flow_nominal=6.36
          "Nominal mass flow rate";

        parameter Modelica.Units.SI.PressureDifference dp1_nominal = 20000
          "Pressure drop at m_flow_nominal";

         parameter Modelica.Units.SI.PressureDifference dp2_nominal = 39989.6
          "Pressure drop at m_flow_nominal";

        Fluid.Sources.Boundary_pT exh(
          redeclare package Medium = MediumAir,
          T(displayUnit="K"),
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{72,-22},{52,-2}})));

        Fluid.Sensors.TemperatureTwoPort Tdra(
          redeclare package Medium = MediumAir,
          m_flow_nominal=m1_flow_nominal,
          tau=30,
          T_start(displayUnit="K")) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={30,-12})));
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
          annotation (Placement(transformation(extent={{-10,-28},{10,-8}})));
        Fluid.Sources.MassFlowSource_T      sou(
          redeclare package Medium = MediumWat,
          use_m_flow_in=true,
          m_flow=m2_flow_nominal,
          use_T_in=true,
          T=303.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={60,-38})));
        Fluid.Sources.MassFlowSource_T
                                  pro1(
          redeclare package Medium = MediumAir,
          use_m_flow_in=true,
          m_flow=m1_flow_nominal,
          use_T_in=true,
          T=573.15,
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{-80,-22},{-60,-2}})));
        Fluid.Sources.Boundary_pT           sou1(
          redeclare package Medium = MediumWat,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-70,-40})));
        Fluid.Sensors.TemperatureTwoPort Tfw(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m2_flow_nominal,
          tau=30,
          T_start(displayUnit="K")) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={-38,-40})));
        Modelica.Blocks.Sources.CombiTimeTable MeaData(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Outputs),
          verboseRead=true,
          columns=2:3,
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-180,-160},{-160,-140}})));
        Modelica.Blocks.Sources.CombiTimeTable Tflu(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={3},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-140,-46},{-120,-26}})));
        Modelica.Blocks.Sources.CombiTimeTable mFloFlu(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={4},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{-176,-14},{-156,6}})));
        Modelica.Blocks.Sources.CombiTimeTable TfwSup(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={5},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{120,-52},{100,-32}})));
        Modelica.Blocks.Sources.CombiTimeTable mFloFwsup(
          tableOnFile=true,
          tableName="table",
          fileName=ModelicaServices.ExternalReferences.loadResource(Inputs),
          verboseRead=true,
          columns={6},
          extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
          timeScale=1)
          annotation (Placement(transformation(extent={{146,-90},{126,-70}})));
        Modelica.Blocks.Sources.RealExpression Tfw_s(y=Tfw.T) "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{-180,-126},{-160,-106}})));
        Modelica.Blocks.Sources.RealExpression Tdr_s(y=Tdra.T) "Fuel heat flow rate"
          annotation (Placement(transformation(extent={{-180,-104},{-160,-84}})));
        Modelica.Blocks.Math.Gain             iQfue1(k=1)
          annotation (Placement(transformation(extent={{-128,-14},{-108,6}})));
        Modelica.Blocks.Math.Gain             iQfue2(k=0.9)
          annotation (Placement(transformation(extent={{106,-90},{86,-70}})));
        Modelica.Blocks.Math.Add diffFW(k1=-1)
          annotation (Placement(transformation(extent={{-120,-120},{-100,-140}})));
        Modelica.Blocks.Math.Add diffDfT(k1=-1)
          annotation (Placement(transformation(extent={{-80,-100},{-60,-120}})));
        parameter Real r_nominal=2/3
          "Ratio between air-side and water-side convective heat transfer (hA-value) at nominal condition";
        parameter Modelica.Units.SI.Temperature T_a2_nominal(displayUnit="K") = 250
          "Nominal temperature at port a2";
        parameter Modelica.Units.SI.Temperature T_a1_nominal(displayUnit="K") = 577
          "Nominal temperature at port a1";
        Modelica.Blocks.Math.Gain diff1(k=-1)
          annotation (Placement(transformation(extent={{-36,-120},{-16,-100}})));
        Modelica.Blocks.Math.Gain diff2(k=0.95)
          annotation (Placement(transformation(extent={{-38,-174},{-18,-154}})));
      equation
        connect(eco.port_a2, sou.ports[1]) annotation (Line(
            points={{10,-24},{20,-24},{20,-38},{50,-38}},
            color={0,127,255},
            thickness=0.5));
        connect(eco.port_b1, Tdra.port_b) annotation (Line(
            points={{10,-12},{20,-12}},
            color={0,0,0},
            thickness=0.5));
        connect(Tdra.port_a, exh.ports[1]) annotation (Line(
            points={{40,-12},{52,-12}},
            color={0,0,0},
            thickness=0.5));
        connect(pro1.ports[1],eco. port_a1) annotation (Line(
            points={{-60,-12},{-10,-12}},
            color={0,0,0},
            thickness=0.5));
        connect(sou1.ports[1], Tfw.port_b) annotation (Line(
            points={{-60,-40},{-48,-40}},
            color={0,127,255},
            thickness=0.5));
        connect(Tfw.port_a, eco.port_b2) annotation (Line(
            points={{-28,-40},{-20,-40},{-20,-24},{-10,-24}},
            color={0,127,255},
            thickness=0.5));
        connect(Tflu.y[1], pro1.T_in) annotation (Line(
            points={{-119,-36},{-90,-36},{-90,-8},{-82,-8}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(TfwSup.y[1], sou.T_in) annotation (Line(
            points={{99,-42},{72,-42}},
            color={0,0,127},
            pattern=LinePattern.Dash));
        connect(mFloFlu.y[1], iQfue1.u)
          annotation (Line(points={{-155,-4},{-130,-4}}, color={0,0,127}));
        connect(iQfue1.y, pro1.m_flow_in)
          annotation (Line(points={{-107,-4},{-82,-4}}, color={0,0,127}));
        connect(mFloFwsup.y[1], iQfue2.u)
          annotation (Line(points={{125,-80},{108,-80}}, color={0,0,127}));
        connect(iQfue2.y, sou.m_flow_in)
          annotation (Line(points={{85,-80},{85,-46},{72,-46}}, color={0,0,127}));
        connect(Tfw_s.y, diffFW.u2) annotation (Line(points={{-159,-116},{-130,-116},{
                -130,-124},{-122,-124}}, color={0,0,127}));
        connect(MeaData.y[1], diffFW.u1) annotation (Line(points={{-159,-150},{-138,-150},
                {-138,-136},{-122,-136}}, color={0,0,127}));
        connect(MeaData.y[2], diffDfT.u1) annotation (Line(points={{-159,-150},{-92,-150},
                {-92,-116},{-82,-116}}, color={0,0,127}));
        connect(Tdr_s.y, diffDfT.u2) annotation (Line(points={{-159,-94},{-100,-94},{-100,
                -104},{-82,-104}}, color={0,0,127}));
        connect(diffDfT.y, diff1.u)
          annotation (Line(points={{-59,-110},{-38,-110}}, color={0,0,127}));
        connect(MeaData.y[1], diff2.u) annotation (Line(points={{-159,-150},{-48,-150},
                {-48,-164},{-40,-164}}, color={0,0,127}));
        annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Components/HeatBalanceBoiler.mos"
              "Simulate and plot"),
          experiment(
            StartTime=30499200,
            StopTime=31016700,
            Tolerance=1e-06,
            __Dymola_Algorithm="Dassl"),            Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
                  {160,100}})),                                        Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})));
      end EconomizerTrain;
    end ExhaustEconomizer;

    package Pumps
      model FWpumpWithPressurisedTanks
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

        Fluid.Storage.ExpansionVessel           tanFW(
          redeclare final package Medium = Medium,
          final V_start=VTanFW_start,
          final p_start=pTanFW)
          "Feedwater tank"
          annotation (Placement(transformation(extent={{-30,-82},{-10,-62}})));
        Fluid.Movers.FlowControlled_m_flow fwPum(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m_flow_nominal,
          addPowerToMedium=false,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal=dp_nominal) "Feed water pump"
          annotation (Placement(transformation(extent={{2,-110},{22,-90}})));
        Fluid.Movers.FlowControlled_m_flow CDpum(
          redeclare package Medium = MediumWat,
          m_flow_nominal=m_flow_nominal,
          addPowerToMedium=false,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal=dp_nominal) "Feed water pump"
          annotation (Placement(transformation(extent={{-80,-110},{-60,-90}})));
        Fluid.Storage.ExpansionVessel           tanFW1(
          redeclare final package Medium = Medium,
          final V_start=VTanFW_start,
          final p_start=pTanFW)
          "Feedwater tank"
          annotation (Placement(transformation(extent={{-110,-82},{-90,-62}})));
        Fluid.Sources.Boundary_pT           sou2(
          redeclare package Medium = MediumWat,
          p=100000,
          T=303.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-158,-100})));
        Fluid.Sources.Boundary_pT           sou3(
          redeclare package Medium = MediumWat,
          p=100000,
          T=303.15,
          nPorts=1)
          "Source"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={90,-100})));
        Fluid.FixedResistances.PressureDrop res(
          redeclare package Medium = MediumAir,
          m_flow_nominal=10,
          dp_nominal=50000)
          annotation (Placement(transformation(extent={{-50,-110},{-30,-90}})));
        Fluid.FixedResistances.PressureDrop res1(
          redeclare package Medium = MediumAir,
          m_flow_nominal=10,
          dp_nominal=50000)
          annotation (Placement(transformation(extent={{-132,-110},{-112,-90}})));
        Fluid.FixedResistances.PressureDrop res2(
          redeclare package Medium = MediumAir,
          m_flow_nominal=10,
          dp_nominal=50000)
          annotation (Placement(transformation(extent={{40,-110},{60,-90}})));
      equation
        connect(tanFW.port_a, fwPum.port_a) annotation (Line(points={{-20,-82},{
                -20,-100},{2,-100}},   color={0,127,255},
            thickness=0.5));
        connect(tanFW1.port_a, CDpum.port_a) annotation (Line(
            points={{-100,-82},{-100,-100},{-80,-100}},
            color={0,127,255},
            thickness=0.5));
        connect(CDpum.port_b, res.port_a) annotation (Line(
            points={{-60,-100},{-50,-100}},
            color={0,127,255},
            thickness=0.5));
        connect(res.port_b, tanFW.port_a) annotation (Line(
            points={{-30,-100},{-20,-100},{-20,-82}},
            color={0,127,255},
            thickness=0.5));
        connect(sou2.ports[1], res1.port_a) annotation (Line(
            points={{-148,-100},{-132,-100}},
            color={0,127,255},
            thickness=0.5));
        connect(res1.port_b, tanFW1.port_a) annotation (Line(
            points={{-112,-100},{-100,-100},{-100,-82}},
            color={0,127,255},
            thickness=0.5));
        connect(fwPum.port_b, res2.port_a) annotation (Line(
            points={{22,-100},{40,-100}},
            color={0,127,255},
            thickness=0.5));
        connect(res2.port_b, sou3.ports[1]) annotation (Line(
            points={{60,-100},{80,-100}},
            color={0,127,255},
            thickness=0.5));
        annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/GEDCalibration/CUBoulder/Components/HeatBalanceBoiler.mos"
              "Simulate and plot"),
          experiment(Tolerance=1e-6, StopTime=3600),Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},
                  {160,100}})),                                        Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{160,100}})));
      end FWpumpWithPressurisedTanks;
    end Pumps;
  end Validation;

  package Examples
    extends Modelica.Icons.ExamplesPackage;

  end Examples;

  package Obsolete
    model BoilerPolynomialExhaustCombustiongas "A equilibrium boiler with water phase change from liquid to vapor, discharging
  saturated steam vapor, with the efficiency curve described by a polynomial."
      extends Buildings.BaseClasses.BaseIconLow;
      extends
        Buildings.Experimental.DHC.BaseClasses.Steam.PartialFourPortInterfaceThreeMedium(
        redeclare final package Medium1 = MediumAir,
        redeclare final package Medium2 = MediumWat,
        redeclare final package Medium3 = MediumSte,
        final T_start=MediumSte.saturationTemperature(p_start));

      // Medium declarations
      replaceable package MediumAir =
          Buildings.Media.CombustionAir
        "Air medium on the combustion side";
      replaceable package MediumWat =
          Buildings.Media.Specialized.Water.TemperatureDependentDensity
        "Water medium - port_a2 (inlet)";
      replaceable package MediumSte = Buildings.Media.Steam
        "Steam medium - port_b2 (oulet)";
      // Initialization
      parameter Boolean fixed_p_start=false "Set to true if p_start is to be used as an explicit initial equation, 
    not an initial guess"     annotation (Dialog(tab="Initialization"));
      // Nominal conditions
      parameter Modelica.Units.SI.PressureDifference dp_nominal(displayUnit="Pa")
        "Pressure drop at nominal mass flow rate"
        annotation (Dialog(group="Nominal condition"));
      parameter Modelica.Units.SI.Power Q_flow_nominal "Nominal heating power";
      parameter Modelica.Units.SI.Temperature T_nominal=373.15 "Temperature used to compute nominal efficiency 
    (only used if efficiency curve depends on temperature)";

      // Efficiency, fuel, and boiler properties
      parameter Buildings.Fluid.Types.EfficiencyCurves effCur=Buildings.Fluid.Types.EfficiencyCurves.Constant
        "Curve used to compute the efficiency";
      parameter Real a[:]={0.9} "Coefficients for efficiency curve";
      parameter Buildings.Fluid.Data.Fuels.Generic fue "Fuel type"
        annotation (choicesAllMatching=true);
      parameter Modelica.Units.SI.ThermalConductance UA=0.05*Q_flow_nominal/30
        "Overall UA value";
      parameter Modelica.Units.SI.Volume V=1.5E-6*Q_flow_nominal
        "Total internal volume of boiler" annotation (Dialog(tab="Dynamics", enable=
             not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

      parameter Modelica.Units.SI.Volume V_com=m1_flow_nominal*tau1/rho1_nominal
        "Total internal volume of combustion" annotation (Dialog(tab="Dynamics",
            enable=not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

      parameter Modelica.Units.SI.Mass mDry=1.5E-3*Q_flow_nominal
        "Mass of boiler that will be lumped to water heat capacity" annotation (
          Dialog(tab="Dynamics", enable=not (energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState)));

      // Variables
      Modelica.Units.SI.Efficiency eta=if effCur == Buildings.Fluid.Types.EfficiencyCurves.Constant
           then a[1] elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.Polynomial
           then Buildings.Utilities.Math.Functions.polynomial(a=a, x=y_internal)
           elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear
           then Buildings.Utilities.Math.Functions.quadraticLinear(
          a=aQuaLin,
          x1=y_internal,
          x2=MediumSte.saturationTemperature(port_a2.p)) else 0 "Boiler efficiency";
      Modelica.Units.SI.Power QFue_flow=y_internal*Q_flow_nominal/eta_nominal
        "Heat released by fuel";

      Modelica.Units.SI.Power QWat_flow=eta*QFue_flow
        "Heat transfer from gas into water";
      Modelica.Units.SI.MassFlowRate mFue_flow=QFue_flow/fue.h
        "Fuel mass flow rate";
      Modelica.Units.SI.VolumeFlowRate VFue_flow=mFue_flow/fue.d
        "Fuel volume flow rate";

      Modelica.Blocks.Interfaces.RealInput y(min=0, max=1)
        "Part load ratio"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealOutput VLiq(
        final quantity="Volume",
        final unit="m3",
        min=0) "Output liquid water volume"
        annotation (Placement(transformation(extent={{100,-90},{120,-70}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
                              "Heat port, can be used to connect to ambient"
        annotation (Placement(transformation(extent={{-10,90},{10,110}})));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heaCapDry(C=500*mDry,
          T(start=T_start)) if not steadyDynamics "Heat capacity of boiler metal"
        annotation (Placement(transformation(extent={{30,-78},{50,-58}})));
      Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation
        vol(
        redeclare final package MediumSte = MediumSte,
        redeclare final package MediumWat = MediumWat,
        final allowFlowReversal=allowFlowReversal2,
        final energyDynamics=energyDynamics,
        final massDynamics=massDynamics,
        final p_start=p_start,
        final m_flow_nominal=m2_flow_nominal,
        final show_T=show_T,
        final V=V,
        final fixed_p_start=fixed_p_start) "Steam/water control volume"
        annotation (Placement(transformation(extent={{8,-50},{28,-30}})));
      Buildings.Fluid.FixedResistances.PressureDrop res(
        redeclare final package Medium = MediumWat,
        final allowFlowReversal=allowFlowReversal2,
        final m_flow_nominal=m2_flow_nominal,
        final show_T=show_T,
        final dp_nominal=dp_nominal) "Flow resistance"
        annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));

      Modelica.Blocks.Interfaces.RealOutput QFueFlo(
        final quantity="HeatFlowRate",
        final unit="W",
        min=0) "Heat flow rate of the fuel"
        annotation (Placement(transformation(extent={{100,80},{120,100}})));

      Fluid.MixingVolumes.MixingVolume fur(
        redeclare package Medium = Medium1,
        m_flow_nominal=m1_flow_nominal,
        V=V_com,
        nPorts=2) "Furnace"
        annotation (Placement(transformation(extent={{6,40},{-14,60}})));
      Fluid.Movers.FlowControlled_m_flow fan(
        redeclare package Medium = Medium1,
        m_flow_nominal=m1_flow_nominal,
        addPowerToMedium=false,
        nominalValuesDefineDefaultPressureCurve=true,
        dp_nominal=10000)
        annotation (Placement(transformation(extent={{-80,50},{-60,30}})));

      Modelica.Blocks.Interfaces.RealOutput mFueFlo(
        final quantity="MassFlowRate",
        final unit="kg/s",
        min=0) "Mass flow rate of the fuel"
        annotation (Placement(transformation(extent={{100,66},{120,86}})));
    protected
      final parameter Boolean steadyDynamics=if energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState
           then true else false "= true, if steady state formulation";
      parameter Real eta_nominal(fixed=false)
        "Boiler efficiency at nominal condition";
      parameter Real aQuaLin[6]=if size(a, 1) == 6 then a else fill(0, 6)
        "Auxiliary variable for efficiency curve because quadraticLinear requires exactly 6 elements";

      Modelica.Blocks.Interfaces.RealInput y_internal(min=0, max=1)
        "Internal block needed for conditional input part load ratio";

      Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo
        if not steadyDynamics "Prescribed heat flow (if heatPort is connected)"
        annotation (Placement(transformation(extent={{-59,-90},{-39,-70}})));
    public
      Modelica.Blocks.Sources.RealExpression Q_flow_in(y=(eta)*port_a1.m_flow*(
            senEntFlo1.H_flow - senEntFlo.H_flow)) if not steadyDynamics
        "Heat transfer from gas into water (if heatPort is connected)"
        annotation (Placement(transformation(extent={{-98,-90},{-78,-70}})));

    protected
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor UAOve(G=UA)
        if not steadyDynamics
        "Overall thermal conductance (if heatPort is connected)"
        annotation (Placement(transformation(extent={{-18,-74},{2,-54}})));

      Modelica.Blocks.Sources.RealExpression QFue_flow_out(y=QFue_flow)
        "Heat flow rate of the fuel"
        annotation (Placement(transformation(extent={{60,80},{80,100}})));
    public
      Modelica.Blocks.Sources.RealExpression Q_flow_exh(y=QFue_flow)
                              "Heat losses in the exhaust"
        annotation (Placement(transformation(extent={{80,50},{60,70}})));
    protected
      HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo1 if not steadyDynamics
        "Prescribed heat flow (if heatPort is connected)"
        annotation (Placement(transformation(extent={{41,50},{21,70}})));
    public
      Modelica.Blocks.Sources.RealExpression m_exh_flow(y=y*((mFueFlo*19.6) +
            mFueFlo))         "Heat losses in the exhaust"
        annotation (Placement(transformation(extent={{-100,-2},{-80,18}})));
      Fluid.Sensors.EnthalpyFlowRate senEntFlo(redeclare package Medium = Medium1,
          m_flow_nominal=m1_flow_nominal)
        annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
      Fluid.Sensors.EnthalpyFlowRate senEntFlo1(redeclare package Medium = Medium1,
          m_flow_nominal=m1_flow_nominal)
        annotation (Placement(transformation(extent={{22,30},{42,50}})));
      Fluid.MixingVolumes.MixingVolume vol2(
        redeclare package Medium = Medium1,
        m_flow_nominal=m1_flow_nominal,
        V=V_com,
        nPorts=2) annotation (Placement(transformation(extent={{60,40},{80,20}})));
    public
      Modelica.Blocks.Sources.RealExpression Q_flow_exh1(y=(1 - eta)*port_a1.m_flow*
            (senEntFlo1.H_flow - senEntFlo.H_flow))
                              "Heat losses in the exhaust"
        annotation (Placement(transformation(extent={{-6,-10},{14,10}})));
    protected
      Modelica.Blocks.Sources.RealExpression mFue_flow_out(y=mFue_flow)
        "Mass flow rate of the fuel"
        annotation (Placement(transformation(extent={{60,66},{80,86}})));
    protected
      HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo2 if not steadyDynamics
        "Prescribed heat flow (if heatPort is connected)"
        annotation (Placement(transformation(extent={{23,-10},{43,10}})));
    initial equation
      if effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear then
        assert(size(a, 1) == 6, "The boiler has the efficiency curve set to 'Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear',
    and hence the parameter 'a' must have exactly 6 elements.
    However, only "     + String(size(a, 1)) + " elements were provided.");
      end if;

      if effCur == Buildings.Fluid.Types.EfficiencyCurves.Constant then
        eta_nominal = a[1];
      elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.Polynomial then
        eta_nominal = Buildings.Utilities.Math.Functions.polynomial(a=a, x=1);
      elseif effCur == Buildings.Fluid.Types.EfficiencyCurves.QuadraticLinear then
        // For this efficiency curve, a must have 6 elements.
        eta_nominal = Buildings.Utilities.Math.Functions.quadraticLinear(
          a=aQuaLin,
          x1=1,
          x2=T_nominal);
      else
        eta_nominal = 999;
      end if;

    equation
      assert(eta > 0.001, "Efficiency curve is wrong.");

      connect(y, y_internal);

      if steadyDynamics then
        -QWat_flow = port_a2.m_flow*actualStream(port_a2.h_outflow) + port_b2.m_flow
          *actualStream(port_b2.h_outflow);
      end if;

      connect(UAOve.port_a, heatPort) annotation (Line(
          points={{-18,-64},{-20,-64},{-20,100},{0,100}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(Q_flow_in.y, preHeaFlo.Q_flow) annotation (Line(
          points={{-77,-80},{-59,-80}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(heaCapDry.port, UAOve.port_b) annotation (Line(points={{40,-78},{40,-80},
              {2,-80},{2,-64}}, color={191,0,0}));
      connect(preHeaFlo.port, UAOve.port_b)
        annotation (Line(points={{-39,-80},{2,-80},{2,-64}}, color={191,0,0}));
      connect(vol.heatPort, UAOve.port_b)
        annotation (Line(points={{18,-50},{18,-64},{2,-64}}, color={191,0,0}));
      connect(vol.port_b, port_b2)
        annotation (Line(points={{28,-40},{100,-40}}, color={0,127,255}));
      connect(port_a2, res.port_a)
        annotation (Line(points={{-100,-40},{-60,-40}}, color={0,127,255}));
      connect(res.port_b, vol.port_a)
        annotation (Line(points={{-40,-40},{8,-40}}, color={0,127,255}));
      connect(QFue_flow_out.y, QFueFlo)
        annotation (Line(points={{81,90},{110,90}}, color={0,0,127}));
      connect(Q_flow_exh.y, preHeaFlo1.Q_flow) annotation (Line(
          points={{59,60},{41,60}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(preHeaFlo1.port, fur.heatPort) annotation (Line(points={{21,60},{18,60},
              {18,50},{6,50}}, color={191,0,0}));
      connect(port_a1, fan.port_a)
        annotation (Line(points={{-100,40},{-80,40}}, color={0,127,255}));
      connect(vol.VLiq, VLiq) annotation (Line(points={{29,-33},{60,-33},{60,-80},{110,
              -80}}, color={0,0,127}));
      connect(m_exh_flow.y, fan.m_flow_in)
        annotation (Line(points={{-79,8},{-70,8},{-70,28}}, color={0,0,127}));
      connect(mFue_flow_out.y, mFueFlo)
        annotation (Line(points={{81,76},{110,76}}, color={0,0,127}));
      connect(fan.port_b, senEntFlo.port_a)
        annotation (Line(points={{-60,40},{-50,40}}, color={0,127,255}));
      connect(senEntFlo.port_b, fur.ports[1])
        annotation (Line(points={{-30,40},{-3,40}}, color={0,127,255}));
      connect(fur.ports[2], senEntFlo1.port_a)
        annotation (Line(points={{-5,40},{22,40}}, color={0,127,255}));
      connect(senEntFlo1.port_b, vol2.ports[1])
        annotation (Line(points={{42,40},{69,40}}, color={0,127,255}));
      connect(vol2.ports[2], port_b1)
        annotation (Line(points={{71,40},{100,40}}, color={0,127,255}));
      connect(Q_flow_exh1.y, preHeaFlo2.Q_flow) annotation (Line(
          points={{15,0},{23,0}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(preHeaFlo2.port, vol2.heatPort)
        annotation (Line(points={{43,0},{54,0},{54,30},{60,30}}, color={191,0,0}));
      annotation (
        defaultComponentName="boi",
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-80,60},{80,-60}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-40,40},{40,-40}},
              fillColor={127,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{20,18},{0,8},{20,-12},{0,-22}},
              color={0,0,0},
              smooth=Smooth.Bezier,
              extent={{-60,-22},{-36,2}}),
            Line(
              points={{-2,18},{-22,8},{-2,-12},{-22,-22}},
              color={0,0,0},
              smooth=Smooth.Bezier,
              extent={{-60,-22},{-36,2}})}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>
This model represents a steam boiler that discharges saturated 
steam and has an efficiency curve defined by a polynomial.
The efficiency in this model represents the fuel-to-water 
efficiency (e.g., thermal efficiency).
This model is similar to the 
<a href=\"modelica://Buildings.Fluid.Boilers.BoilerPolynomial\"> 
Buildings.Fluid.Boilers.BoilerPolynomial</a> for the efficiency 
and fuel mass flow rate computation with the following exceptions:
</p>
<ul>
<li>
Water enters <code>port_a</code> in liquid state and exits 
<code>port_b</code> in vapor state.
</li> 
<li>
The liquid and vapor phases are at equilibrium; thus, the steam
boiler is constrained to saturated states only with the volume
containing a wet steam mixture. 
</li>
<li>
If the boiler is configured in steady state, several blocks involving
the heat flow rate are conditionally removed to avoid overconstraining
the model. This is because the discharging fluid is constrained at 
a saturated state. The blocks that are conditionally removed as a 
result are within the green region in the below figure:
</li>
</ul>

<p align=\"center\">
<img src=\"modelica://Buildings/Resources/Images/Experimental/DHC/Plants/Steam/BaseClasses/BoilerPolynomial.png\" border=\"1\"
alt=\"Boiler polynomial steam with blocks in green conditionally removed if steady state\"/>
</p>
<h4>Implementation</h4>
<p>
In order to improve the numerical efficiency, this model follows 
the split-medium approach using the
<a href=\"modelica://Buildings.Fluid.Interfaces.PartialTwoPortTwoMedium\">
Buildings.Fluid.Interfaces.PartialTwoPortTwoMedium</a> interface model.
The saturated mixing volume for an evaporation process 
<a href=\"modelica://Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation\">
Buildings.Experimental.DHC.Plants.Steam.BaseClasses.ControlVolumeEvaporation</a> 
represents the phase change process of water from liquid 
to vapor at equilibrium.
</p>
<h4>Reference</h4>
<p>
Hinkelman, Kathryn, Saranya Anbarasu, Michael Wetter, 
Antoine Gautier, and Wangda Zuo. 2022. “A Fast and Accurate Modeling 
Approach for Water and Steam Thermodynamics with Practical 
Applications in District Heating System Simulation.” Preprint. February 24. 
<a href=\"http://dx.doi.org/10.13140/RG.2.2.20710.29762\">doi:10.13140/RG.2.2.20710.29762</a>.
</p>
</html>",     revisions="<html>
<ul>
<li>
February 25, 2022 by Kathryn Hinkelman:<br/>
Refactored base classes for improved extensibility and relocated models into Steam subpackages.
</li>
<li>
July 22, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"));
    end BoilerPolynomialExhaustCombustiongas;
  end Obsolete;

  package BaseClasses
    block FurnaceHeatBalance

    parameter Modelica.Units.SI.Temperature T_exh_nominal=373.15 "Exhaust temperature used to compute nominal efficiency";
      Modelica.Blocks.Interfaces.RealOutput eta_boi "Efficiency of boiler"
        annotation (Placement(transformation(extent={{100,10},{120,30}}),
            iconTransformation(extent={{100,10},{120,30}})));
      Modelica.Blocks.Interfaces.RealOutput per_exh "Exhaust heat percentage"
        annotation (Placement(transformation(extent={{100,-30},{120,-10}}),
            iconTransformation(extent={{100,-40},{120,-20}})));
      Modelica.Blocks.Tables.CombiTable1Ds hCO2(
        table=[200,-3.41; 293,-0.19; 298,0.0; 400,4.01; 500,8.31; 600,47.32; 700,17.76],
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
        "Enthalpy of carbondioxide at the T exhaust (h in MJ/kg mol)"
        annotation (Placement(transformation(extent={{-20,60},{0,80}})));

      Modelica.Blocks.Tables.CombiTable1Ds hO2(
        table=[200,-2.87; 293,-0.15; 298,0; 400,3.03; 500,6.09; 600,9.25; 700,12.50],
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
        "Enthalpy of oxygen at the T exhaust "
        annotation (Placement(transformation(extent={{-20,30},{0,50}})));

      Modelica.Blocks.Tables.CombiTable1Ds hN2(
        table=[200,-2.86; 293,-0.15; 298,0; 400,2.97; 500,5.91; 600,8.89; 700,11.94],
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
                                       "Enthalpy of nitrogen at the T exhaust "
        annotation (Placement(transformation(extent={{-20,0},{0,20}})));

      Modelica.Blocks.Sources.RealExpression Texh(y=T_exh_nominal)
        "Exhaust temperature setpoint"
        annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
      Modelica.Blocks.Interfaces.RealInput qLos "Heat loss from the boiler casing"
        annotation (Placement(transformation(extent={{-140,40},{-100,80}}),
            iconTransformation(extent={{-140,40},{-100,80}})));
      parameter Modelica.Units.SI.SpecificEnthalpy QFue "Heating value of per kg mol of fuel";
      parameter Real FA_ratio "Fuel air ratio, alpha (20% excess air, FA_ratio = 1.20)";

      //Variables
      Modelica.Units.SI.SpecificEnthalpy QUse "Useful heat from the combustion process";
      Modelica.Units.SI.MolarMass nCO2 = 1 "Moles of CO2";
      Modelica.Units.SI.MolarMass nH2O = 2 "Moles of H2O";
      Modelica.Units.SI.MolarMass nO2 = 2*(FA_ratio-1) "Moles of O2";
      Modelica.Units.SI.MolarMass nN2 = 2* FA_ratio * 3.76 "Moles of N2";
      Modelica.Units.SI.MolarMass nCH4 = 1 "Moles of Fuel, CH4";
      Modelica.Units.SI.MolarMass nO2N2 = 2*FA_ratio "Moles of air";

      //Enthalpy
     Modelica.Units.SI.SpecificEnthalpy h_CO2 = ((nCO2/nCH4)* hCO2.y[1]*10^6);
     Modelica.Units.SI.SpecificEnthalpy h_H2O = ((nH2O/nCH4)* hH2O.y[1]*10^6);
     Modelica.Units.SI.SpecificEnthalpy h_O2 = ((nO2/nCH4)* hO2.y[1]*10^6);
     Modelica.Units.SI.SpecificEnthalpy h_N2 = ((nN2/nCH4)* hN2.y[1]*10^6);
     Modelica.Units.SI.SpecificEnthalpy ha = ((nO2N2/nCH4)* (hN1.y[1]+hO1.y[1])*10^6);

     Real fs = (nCH4*16)/((29*(nCO2+ (nH2O*2)/4 - (2*nO2)/2))*4.76)  "Fuel air ratio by weight";

     //test outputs
     Real per_exh_prod;
     Real per_wat;
     Real per_inletair;
     Real per_total;
     Real per_loss;

      Modelica.Blocks.Tables.CombiTable1Ds hH2O(
        table=[200,-3.28; 293,-0.17; 298,0; 400,3.45; 500,6.92; 600,10.50; 700,14.18],
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
        "Enthalpy of Water vapor at the T exhaust "
        annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));

      Modelica.Blocks.Tables.CombiTable1Ds hO1(table=[200,-2.87; 293,-0.15; 298,0; 400,
            3.03; 500,6.09; 600,9.25],
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
        "Enthalpy of oxygen at the T exhaust "
        annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));
      Modelica.Blocks.Tables.CombiTable1Ds hN1(table=[200,-2.86; 293,-0.15; 298,0; 400,
            2.97; 500,5.91; 600,8.89],
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
                                       "Enthalpy of nitrogen at the T exhaust "
        annotation (Placement(transformation(extent={{-20,-90},{0,-70}})));
      Modelica.Blocks.Interfaces.RealInput T_Air_in "Enthalpy of air" annotation (
          Placement(transformation(extent={{-140,-80},{-100,-40}}),
            iconTransformation(extent={{-140,-40},{-100,0}})));
    equation

      QUse = QFue*16 + ha- h_CO2-h_H2O-h_O2-h_N2-(nH2O*2260000*18);
      eta_boi = QUse/(QFue*16);
      per_exh = ((QFue*16)-QUse-(nH2O*2260000*18))/(QFue*16);
      //per_exh = ((QFue*16)-QUse)/(QFue*16);

      per_exh_prod= (h_CO2+h_H2O+h_O2+h_N2)/(QFue*16);
      per_wat = (nH2O*2260000*18)/(QFue*16);
      per_inletair= ha/(QFue*16);
      per_loss = qLos/(QFue*16);
      per_total = eta_boi+per_exh+per_wat;

      connect(Texh.y, hCO2.u)
        annotation (Line(points={{-59,70},{-22,70}}, color={0,0,127},
          pattern=LinePattern.Dash));
      connect(Texh.y, hO2.u) annotation (Line(points={{-59,70},{-48,70},{-48,40},
              {-22,40}},
                    color={0,0,127},
          pattern=LinePattern.Dash));
      connect(Texh.y, hN2.u) annotation (Line(points={{-59,70},{-48,70},{-48,10},
              {-22,10}},
                    color={0,0,127},
          pattern=LinePattern.Dash));
      connect(Texh.y, hH2O.u) annotation (Line(points={{-59,70},{-48,70},{-48,-20},
              {-22,-20}},color={0,0,127},
          pattern=LinePattern.Dash));
      connect(T_Air_in, hO1.u) annotation (Line(points={{-120,-60},{-46,-60},{-46,
              -50},{-22,-50}},
                          color={0,0,127},
          pattern=LinePattern.Dash));
      connect(T_Air_in, hN1.u) annotation (Line(points={{-120,-60},{-46,-60},{-46,
              -80},{-22,-80}},
                          color={0,0,127},
          pattern=LinePattern.Dash));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end FurnaceHeatBalance;

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
  end BaseClasses;
end Components;
