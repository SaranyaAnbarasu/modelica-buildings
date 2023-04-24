within Buildings.GEDCalibration.CUBoulder;
package Components

  block FurnaceHeatBalance

  parameter Modelica.Units.SI.Temperature T_exh_nominal=373.15 "Exhaust temperature used to compute nominal efficiency";
    Modelica.Blocks.Interfaces.RealOutput eta_boi "Efficiency of boiler"
      annotation (Placement(transformation(extent={{100,10},{120,30}}),
          iconTransformation(extent={{100,10},{120,30}})));
    Modelica.Blocks.Interfaces.RealOutput per_exh "Exhaust heat percentage"
      annotation (Placement(transformation(extent={{100,-30},{120,-10}}),
          iconTransformation(extent={{100,-40},{120,-20}})));
    Modelica.Blocks.Tables.CombiTable1Ds hCO2(table=[200,-3.41; 293,-0.19; 298,0.0;
          400,4.01; 500,8.31; 600,47.32])
      "Enthalpy of carbondioxide at the T exhaust (h in MJ/kg mol)"
      annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
    Modelica.Blocks.Tables.CombiTable1Ds hO2(table=[200,-2.87; 293,-0.15; 298,0; 400,
          3.03; 500,6.09; 600,9.25])
      "Enthalpy of oxygen at the T exhaust "
      annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
    Modelica.Blocks.Tables.CombiTable1Ds hN2(table=[200,-2.86; 293,-0.15; 298,0; 400,
          2.97; 500,5.91; 600,8.89]) "Enthalpy of nitrogen at the T exhaust "
      annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
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

    Modelica.Blocks.Tables.CombiTable1Ds hH2O(table=[200,-3.28; 293,-0.17; 298,0;
          400,3.45; 500,6.92; 600,10.50])
      "Enthalpy of Water vapor at the T exhaust "
      annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
    Modelica.Blocks.Tables.CombiTable1Ds hO1(table=[200,-2.87; 293,-0.15; 298,0; 400,
          3.03; 500,6.09; 600,9.25])
      "Enthalpy of oxygen at the T exhaust "
      annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
    Modelica.Blocks.Tables.CombiTable1Ds hN1(table=[200,-2.86; 293,-0.15; 298,0; 400,
          2.97; 500,5.91; 600,8.89]) "Enthalpy of nitrogen at the T exhaust "
      annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
    Modelica.Blocks.Interfaces.RealInput T_Air_in "Enthalpy of air" annotation (
        Placement(transformation(extent={{-140,-80},{-100,-40}}),
          iconTransformation(extent={{-140,-40},{-100,0}})));
  equation

    QUse = QFue*16 + ha- h_CO2-h_H2O-h_O2-h_N2-(nH2O*2260000*18)-qLos;
    eta_boi = QUse/(QFue*16);
    per_exh = ((QFue*16)-QUse-qLos)/(QFue*16);

    connect(Texh.y, hCO2.u)
      annotation (Line(points={{-59,70},{-42,70}}, color={0,0,127}));
    connect(Texh.y, hO2.u) annotation (Line(points={{-59,70},{-48,70},{-48,40},{-42,
            40}}, color={0,0,127}));
    connect(Texh.y, hN2.u) annotation (Line(points={{-59,70},{-48,70},{-48,10},{-42,
            10}}, color={0,0,127}));
    connect(Texh.y, hH2O.u) annotation (Line(points={{-59,70},{-48,70},{-48,-20},{
            -42,-20}}, color={0,0,127}));
    connect(T_Air_in, hO1.u) annotation (Line(points={{-120,-60},{-50,-60},{-50,-50},
            {-42,-50}}, color={0,0,127}));
    connect(T_Air_in, hN1.u) annotation (Line(points={{-120,-60},{-50,-60},{-50,-80},
            {-42,-80}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end FurnaceHeatBalance;

  model BoilerPolynomialExhaust "A equilibrium boiler with water phase change from liquid to vapor, discharging
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
  end BoilerPolynomialExhaustCombustiongas;

  model BoilerPolynomialExhaustHeatBalance "A equilibrium boiler with water phase change from liquid to vapor, discharging
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
    Modelica.Units.SI.Power QWat_flow=furnaceHeatBalance.eta_boi*QFue_flow
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
      redeclare package Medium = MediumAir,
      energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
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
    Modelica.Blocks.Sources.RealExpression m_exh_flow(y=y*((mFueFlo/
          furnaceHeatBalance.fs) + mFueFlo))
                    "Heat losses in the exhaust"
      annotation (Placement(transformation(extent={{-100,-2},{-80,18}})));
    Fluid.Sensors.TemperatureTwoPort
                                   senTem(   redeclare package Medium = MediumAir,
        m_flow_nominal=m1_flow_nominal,
      T_start(displayUnit="K"))
      annotation (Placement(transformation(extent={{-48,30},{-28,50}})));
    FurnaceHeatBalance furnaceHeatBalance(
      T_exh_nominal=T_exh_nominal,        QFue=fue.h, FA_ratio=FA_ratio)
      annotation (Placement(transformation(extent={{-60,80},{-40,100}})));
  public
    Modelica.Blocks.Sources.RealExpression qLos(y=heatPort.Q_flow)
      "Heat losses ifrom boiler casing"
      annotation (Placement(transformation(extent={{-102,86},{-82,106}})));
    parameter Real FA_ratio
      "Fuel air ratio, alpha (20% excess air, FA_ratio = 1.20)";
    parameter Modelica.Units.SI.Temperature T_exh_nominal=373.15
      "Exhaust temperature used to compute nominal efficiency";
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
      annotation (Line(points={{-100,40},{-80,40}}, color={0,127,255}));
    connect(vol.VLiq, VLiq) annotation (Line(points={{29,-33},{60,-33},{60,-80},{110,
            -80}}, color={0,0,127}));
    connect(m_exh_flow.y, fan.m_flow_in)
      annotation (Line(points={{-79,8},{-70,8},{-70,28}}, color={0,0,127}));
    connect(mFue_flow_out.y, mFueFlo)
      annotation (Line(points={{81,76},{110,76}}, color={0,0,127}));
    connect(fan.port_b, senTem.port_a)
      annotation (Line(points={{-60,40},{-48,40}}, color={0,127,255}));
    connect(senTem.port_b, vol1.ports[2])
      annotation (Line(points={{-28,40},{1,40}}, color={0,127,255}));
    connect(senTem.T, furnaceHeatBalance.T_Air_in) annotation (Line(points={{-38,51},
            {-38,70},{-80,70},{-80,88},{-62,88}}, color={0,0,127}));
    connect(qLos.y, furnaceHeatBalance.qLos)
      annotation (Line(points={{-81,96},{-62,96}}, color={0,0,127}));
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
  end BoilerPolynomialExhaustHeatBalance;

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

  package Validation
    extends Modelica.Icons.ExamplesPackage;

    model BoilerPolynomialFourPortHeatBalance
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
      BoilerPolynomialExhaustHeatBalance
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
      Modelica.Blocks.Sources.Pulse pulse(
        amplitude=1,
        width=30,
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
      Modelica.Blocks.Sources.RealExpression QFlue(y=boiDyn.Q_flow_exh.y)
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
      connect(pulse.y, boiDyn.y) annotation (Line(points={{-99,124},{-40,124},{-40,
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
</html>",     revisions="<html>
<ul>
<li>
July 23, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"),
        Diagram(coordinateSystem(extent={{-180,-100},{180,160}}), graphics={
              Rectangle(
              extent={{-180,-40},{180,-180}},
              lineColor={28,108,200},
              fillColor={174,179,179},
              fillPattern=FillPattern.Solid), Text(
              extent={{-176,-26},{-70,-46}},
              textColor={28,108,200},
              textString="Steady state verification (Heat balance)")}));
    end BoilerPolynomialFourPortHeatBalance;

    model FurnaceHeatBalance
        extends Modelica.Icons.Example;
      Buildings.GEDCalibration.CUBoulder.Components.FurnaceHeatBalance
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

    model BoilerPolynomialFourPort
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
      GEDCalibration.CUBoulder.Components.BoilerPolynomialExhaust boiDyn(
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
</html>",     revisions="<html>
<ul>
<li>
July 23, 2021 by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"),
        Diagram(coordinateSystem(extent={{-180,-100},{180,140}}), graphics={
              Rectangle(
              extent={{-180,-40},{180,-180}},
              lineColor={28,108,200},
              fillColor={174,179,179},
              fillPattern=FillPattern.Solid), Text(
              extent={{-176,-26},{-70,-46}},
              textColor={28,108,200},
              textString="Steady state verification (Heat balance)")}));
    end BoilerPolynomialFourPort;

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
        annotation (Placement(transformation(extent={{-64,-54},{-44,-34}})));
      Fluid.Movers.FlowControlled_m_flow fwPum(
        redeclare package Medium = MediumWat,
        m_flow_nominal=m_flow_nominal,
        addPowerToMedium=false,
        nominalValuesDefineDefaultPressureCurve=true,
        dp_nominal=dp_nominal) "Feed water pump"
        annotation (Placement(transformation(extent={{-24,-54},{-4,-34}})));
      Modelica.Blocks.Sources.Constant
                                   const1(k=5.98)
        annotation (Placement(transformation(extent={{-140,-30},{-120,-10}})));
      Modelica.Blocks.Sources.Constant
                                   const2(k=427)
        annotation (Placement(transformation(extent={{-140,-70},{-120,-50}})));
      Modelica.Blocks.Sources.Constant
                                   const3(k=316)
        annotation (Placement(transformation(extent={{-104,10},{-84,30}})));
      Fluid.Sources.Boundary_pT           sou1(
        redeclare package Medium = MediumAir,
        p=110000,
        use_T_in=true,
        T=573.15,
        nPorts=1)
        "Source"
        annotation (Placement(transformation(extent={{-64,6},{-44,26}})));
      HeatTransfer.Sources.FixedTemperature           TAmb(T=288.15)
        "Ambient temperature in boiler room"
        annotation (Placement(transformation(extent={{-14,38},{6,58}})));
      Fluid.Sources.Boundary_pT           sou2(redeclare package Medium = MediumAir,

        p=110000,
        T(displayUnit="K") = 439,
          nPorts=1)
        "Source"
        annotation (Placement(transformation(extent={{96,6},{76,26}})));
      Fluid.Sources.Boundary_pT           sin(
        redeclare package Medium = MediumSte,
        p(displayUnit="bar") = 900000,
        T=453.15,
        nPorts=1)
        "Sink"
        annotation (Placement(transformation(extent={{96,-54},{76,-34}})));
      BoilerPolynomialExhaustHeatBalance boi(
        m1_flow_nominal=7,
        m2_flow_nominal=m_flow_nominal,
        p_start=900000,
        redeclare package MediumAir = MediumAir,
        redeclare package MediumWat = MediumWat,
        redeclare package MediumSte = MediumSte,
        dp_nominal=dp_nominal,
        Q_flow_nominal=Q_flow_nominal,
        fue=Buildings.Fluid.Data.Fuels.Generic(
            h=47522000,
            d=0.8,
            mCO2=2.25),
        UA=0.05*Q_flow_nominal/100,
        V=8.76,
        V_com=10,
        FA_ratio=1.15,
        T_exh_nominal(displayUnit="K") = 439)
        annotation (Placement(transformation(extent={{12,-10},{32,10}})));
      Modelica.Blocks.Sources.Ramp  ramp(duration=3600)
        annotation (Placement(transformation(extent={{-108,50},{-88,70}})));
    equation
      connect(sou.ports[1],fwPum. port_a)
        annotation (Line(points={{-44,-44},{-24,-44}}, color={0,127,255}));
      connect(const1.y,fwPum. m_flow_in)
        annotation (Line(points={{-119,-20},{-14,-20},{-14,-32}},color={0,0,127}));
      connect(const2.y, sou.T_in) annotation (Line(points={{-119,-60},{-74,-60},{-74,
              -40},{-66,-40}}, color={0,0,127}));
      connect(const3.y, sou1.T_in)
        annotation (Line(points={{-83,20},{-66,20}}, color={0,0,127}));
      connect(fwPum.port_b, boi.port_a2) annotation (Line(points={{-4,-44},{6,-44},{
              6,-4},{12,-4}}, color={0,127,255}));
      connect(boi.port_b2, sin.ports[1]) annotation (Line(points={{32,-4},{70,-4},{70,
              -44},{76,-44}}, color={0,127,255}));
      connect(boi.port_b1, sou2.ports[1]) annotation (Line(points={{32,4},{50,4},{50,
              16},{76,16}}, color={0,127,255}));
      connect(sou1.ports[1], boi.port_a1) annotation (Line(points={{-44,16},{-14,16},
              {-14,4},{12,4}}, color={0,127,255}));
      connect(TAmb.port, boi.heatPort) annotation (Line(points={{6,48},{18,48},{18,10},
              {22,10}}, color={191,0,0}));
      connect(ramp.y, boi.y) annotation (Line(points={{-87,60},{-68,60},{-68,66},
              {-20,66},{-20,8},{10,8}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},
                {160,100}})),                                        Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},{160,100}})));
    end HeatBalanceBoiler;
  end Validation;

  package Examples
    extends Modelica.Icons.ExamplesPackage;

  end Examples;

end Components;
