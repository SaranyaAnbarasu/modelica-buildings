within Buildings.Experimental.DHC.BaseClasses.Steam;
partial model PartialFourPortInterfaceThreeMedium
  "Partial model transporting fluid between two ports without storing mass or energy"

//Partial fourport
  replaceable package Medium1 =
    Modelica.Media.Interfaces.PartialMedium "Medium model for combustion side";

  replaceable package Medium2 =
    Modelica.Media.Interfaces.PartialMedium "Medium model for port_a (inlet)";

  replaceable package Medium3 =
    Modelica.Media.Interfaces.PartialMedium "Medium model for port_a (inlet)";

//Fourport two medium
//Massflow rates and small mass flow rates
  parameter Modelica.Units.SI.MassFlowRate m1_flow_nominal(min=0)
    "Nominal mass flow rate" annotation (Dialog(group="Nominal condition"));
  parameter Modelica.Units.SI.MassFlowRate m2_flow_nominal(min=0)
    "Nominal mass flow rate" annotation (Dialog(group="Nominal condition"));
  parameter Modelica.Units.SI.MassFlowRate m1_flow_small(min=0) = 1E-4*abs(m1_flow_nominal)
    "Small mass flow rate for regularization of zero flow"
    annotation(Dialog(tab = "Advanced"));
  parameter Modelica.Units.SI.MassFlowRate m2_flow_small(min=0) = 1E-4*abs(m2_flow_nominal)
    "Small mass flow rate for regularization of zero flow"
    annotation(Dialog(tab = "Advanced"));

  // Diagnostics
  parameter Boolean show_T = false
    "= true, if actual temperature at port is computed"
    annotation (
      Dialog(tab="Advanced", group="Diagnostics"),
      HideResult=true);

 //Starts
  Modelica.Units.SI.MassFlowRate m1_flow = port_a1.m_flow
    "Mass flow rate from port_a1 to port_b1 (m1_flow > 0 is design flow direction)";
  Modelica.Units.SI.PressureDifference dp1(displayUnit="Pa") = port_a1.p - port_b1.p
    "Pressure difference between port_a1 and port_b1";


  Modelica.Units.SI.MassFlowRate m2_flow(start=_m_flow_start) = port_a2.m_flow
    "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";

  Modelica.Units.SI.PressureDifference dp2(start=_dp_start, displayUnit="Pa") = port_a2.p - port_b2.p
    "Pressure difference between port_a and port_b";

//Medium states
  Medium1.ThermodynamicState sta_a1=
    if allowFlowReversal1 then
      Medium1.setState_phX(port_a1.p,
                          noEvent(actualStream(port_a1.h_outflow)),
                          noEvent(actualStream(port_a1.Xi_outflow)))
    else
      Medium1.setState_phX(port_a1.p,
                          inStream(port_a1.h_outflow),
                          inStream(port_a1.Xi_outflow))
      if show_T "Medium properties in port_a1";
  Medium1.ThermodynamicState sta_b1=
    if allowFlowReversal1 then
      Medium1.setState_phX(port_b1.p,
                          noEvent(actualStream(port_b1.h_outflow)),
                          noEvent(actualStream(port_b1.Xi_outflow)))
    else
      Medium1.setState_phX(port_b1.p,
                          port_b1.h_outflow,
                          port_b1.Xi_outflow)
       if show_T "Medium properties in port_b1";


  //Split medium
  Medium2.ThermodynamicState sta_a2=
      Medium2.setState_phX(port_a2.p,
                          noEvent(actualStream(port_a2.h_outflow)),
                          noEvent(actualStream(port_a2.Xi_outflow)))
                                                                    if show_T
                          "Medium properties in port_a2";

  Medium3.ThermodynamicState sta_b2=
      Medium3.setState_phX(port_b2.p,
                          noEvent(actualStream(port_b2.h_outflow)),
                          noEvent(actualStream(port_b2.Xi_outflow)))
                                                                    if show_T
                          "Medium properties in port_b2";


//Defining ports

  Modelica.Fluid.Interfaces.FluidPort_a port_a1(
                     redeclare final package Medium = Medium1,
                     m_flow(min=if allowFlowReversal1 then -Modelica.Constants.inf else 0),
                     h_outflow(start = Medium1.h_default, nominal = Medium1.h_default))
    "Fluid connector a1 (positive design flow direction is from port_a1 to port_b1)"
    annotation (Placement(transformation(extent={{-110,30},{-90,50}})));
  Modelica.Fluid.Interfaces.FluidPort_b port_b1(
                     redeclare final package Medium = Medium1,
                     m_flow(max=if allowFlowReversal1 then +Modelica.Constants.inf else 0),
                     h_outflow(start = Medium1.h_default, nominal = Medium1.h_default))
    "Fluid connector b1 (positive design flow direction is from port_a1 to port_b1)"
    annotation (Placement(transformation(extent={{110,30},{90,50}})));

  Modelica.Fluid.Interfaces.FluidPort_a port_a2(
                     redeclare final package Medium = Medium2,
                     m_flow(min=if allowFlowReversal2 then -Modelica.Constants.inf else 0),
                     h_outflow(start = Medium2.h_default, nominal = Medium2.h_default))
    "Fluid connector a2 (positive design flow direction is from port_a2 to port_b2)"
    annotation (Placement(transformation(extent={{-110,-50},{-90,-30}})));
  Modelica.Fluid.Interfaces.FluidPort_b port_b2(
                     redeclare final package Medium = Medium3,
                     m_flow(max=if allowFlowReversal2 then +Modelica.Constants.inf else 0),
                     h_outflow(start = Medium3.h_default, nominal = Medium3.h_default))
    "Fluid connector b2 (positive design flow direction is from port_a2 to port_b2)"
    annotation (Placement(transformation(extent={{110,-50},{90,-30}})));

  //Medium density



// Flow reversals
  parameter Boolean allowFlowReversal1 = true
    "= false to simplify equations, assuming, but not enforcing, no flow reversal for medium 1"
    annotation(Dialog(tab="Assumptions"), Evaluate=true);
  parameter Boolean allowFlowReversal2 = true
    "= false to simplify equations, assuming, but not enforcing, no flow reversal for medium 2"
    annotation(Dialog(tab="Assumptions"), Evaluate=true);


 //Dynamics
  parameter Modelica.Fluid.Types.Dynamics energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
    "Type of energy balance: dynamic (3 initialization options) or steady state"
    annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));

  parameter Modelica.Fluid.Types.Dynamics massDynamics=energyDynamics
    "Type of mass balance: dynamic (3 initialization options) or steady state"
    annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));

  // Initialization
  parameter Medium3.AbsolutePressure p_start = Medium3.p_default
    "Start value of pressure"
    annotation(Dialog(tab = "Initialization"));
  parameter Medium3.Temperature T_start=Medium3.T_default
    "Start value of temperature"
    annotation(Dialog(tab = "Initialization"));

protected
  final parameter Modelica.Units.SI.MassFlowRate _m_flow_start = 0
  "Start value for m_flow, used to avoid a warning if not set in m_flow, and to avoid m_flow.start in parameter window";
  final parameter Modelica.Units.SI.PressureDifference _dp_start(displayUnit="Pa") = 0
  "Start value for dp, used to avoid a warning if not set in dp, and to avoid dp.start in parameter window";


  parameter Modelica.Units.SI.Time tau1=30 "Time constant at nominal flow"
    annotation (Dialog(tab="Dynamics", group="Nominal condition"));

 parameter Medium1.ThermodynamicState sta1_nominal=Medium1.setState_pTX(
      T=Medium1.T_default, p=Medium1.p_default, X=Medium1.X_default);
  parameter Modelica.Units.SI.Density rho1_nominal=Medium1.density(sta1_nominal)
    "Density, used to compute fluid volume";

  Medium1.ThermodynamicState state_a1_inflow=
  Medium1.setState_phX(port_a1.p, inStream(port_a1.h_outflow), inStream(port_a1.Xi_outflow))
  "state for medium inflowing through port_a1";
  Medium1.ThermodynamicState state_b1_inflow=
  Medium1.setState_phX(port_b1.p, inStream(port_b1.h_outflow), inStream(port_b1.Xi_outflow))
  "state for medium inflowing through port_b1";

  annotation (
  preferredView="info",
    Documentation(info="<html>
<p>
This component defines the interface for models that
transport two fluid streams between four ports.
It is similar to
<a href=\"modelica://Buildings.Fluid.Interfaces.PartialTwoPortInterface\">
Buildings.Fluid.Interfaces.PartialTwoPortInterface</a>,
but it has four ports instead of two.
</p>
<p>
The model is used by other models in this package that add heat transfer,
mass transfer and pressure drop equations.
</p>
</html>", revisions="<html>
<ul>
<li>
February 3, 2022, by Michael Wetter:<br/>
If <code>allowFlowReversal==false</code>, removed <code>noEvent()</code> declaration
for <code>sta_a</code> and for <code>sta_b</code> because the variable is either
already used with <code>inStream()</code> in the computation of <code>state_*_inflow</code>,
or the result of a variable of the model that already may generate an event.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1578\">IBPSA, #1578</a>.
</li>
<li>
February 2, 2022, by Hongxiang Fu:<br/>
If <code>allowFlowReversal==false</code>, replaced <code>actualStream()</code>
with <code>inStream()</code> for <code>sta_a</code> and
removed <code>actualStream()</code> for <code>sta_b</code>.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1578\">IBPSA, #1578</a>.
</li>
<li>
March 30, 2021, by Michael Wetter:<br/>
Added annotation <code>HideResult=true</code>.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1459\">IBPSA, #1459</a>.
</li>
<li>
November 3, 2016, by Michael Wetter:<br/>
Moved computation of pressure drop to variable assignment so that
the model won't mix graphical with textual modeling if used as a base
class for a graphically implemented model.
</li>
<li>
November 3, 2016, by Michael Wetter:<br/>
Removed start values for mass flow rate and pressure difference
to simplify the parameter window.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/552\">#552</a>.
</li>
<li>
January 22, 2016, by Michael Wetter:<br/>
Corrected type declaration of pressure difference.
This is
for <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/404\">#404</a>.
</li>
<li>
November 13, 2013 by Michael Wetter:<br/>
Removed assignment of <code>min</code> and <code>max</code>
attributes of port mass flow rates, as this is already
done in the base class.
</li>
<li>
November 12, 2013 by Michael Wetter:<br/>
Removed <code>import Modelica.Constants;</code> statement.
</li>
<li>
November 11, 2013 by Michael Wetter:<br/>
Removed the parameter <code>homotopyInitialization</code>
as it is no longer used in this model.
</li>
<li>
November 10, 2013 by Michael Wetter:<br/>
In the computation of <code>sta_a1</code>,
<code>sta_a2</code>, <code>sta_b1</code> and <code>sta_b2</code>,
removed the branch that uses the homotopy operator.
The rational is that these variables are conditionally enables (because
of <code>... if show_T</code>). Therefore, the Modelica Language Specification
does not allow for these variables to be used in any equation. Hence,
the use of the homotopy operator is not needed here.
</li>
<li>
October 10, 2013 by Michael Wetter:<br/>
Added <code>noEvent</code> to the computation of the states at the port.
This is correct, because the states are only used for reporting, but not
to compute any other variable.
Use of the states to compute other variables would violate the Modelica
language, as conditionally removed variables must not be used in any equation.
</li>
<li>
October 8, 2013 by Michael Wetter:<br/>
Removed the computation of <code>V_flow</code> and removed the parameter
<code>show_V_flow</code>.
The reason is that the computation of <code>V_flow</code> required
the use of <code>sta_a</code> (to compute the density),
but <code>sta_a</code> is also a variable that is conditionally
enabled. However, this was not correct Modelica syntax as conditional variables
can only be used in a <code>connect</code>
statement, not in an assignment. Dymola 2014 FD01 beta3 is checking
for this incorrect syntax. Hence, <code>V_flow</code> was removed as its
conditional implementation would require a rather cumbersome implementation
that uses a new connector that carries the state of the medium.
</li>
<li>
April 26, 2013 by Marco Bonvini:<br/>
Moved the definitions of <code>dp1</code> and <code>dp2</code> because they cause some problem with PyFMI.
</li>
<li>
March 27, 2012 by Michael Wetter:<br/>
Replaced the erroneous function call <code>Medium.density</code> with
<code>Medium1.density</code> and <code>Medium2.density</code>.
Changed condition to remove <code>sta_a1</code> and <code>sta_a2</code> to also
compute the states at the inlet port if <code>show_V_flow=true</code>.
The previous implementation resulted in a translation error
if <code>show_V_flow=true</code>, but worked correctly otherwise
because the erroneous function call is removed if  <code>show_V_flow=false</code>.
</li>
<li>
March 27, 2011 by Michael Wetter:<br/>
Added <code>homotopy</code> operator.
</li>
<li>
March 21, 2010 by Michael Wetter:<br/>
Changed pressure start value from <code>system.p_start</code>
to <code>Medium.p_default</code> since HVAC models may have water and
air, which are typically at different pressures.
</li>
<li>
September 19, 2008 by Michael Wetter:<br/>
Added equations for the mass balance of extra species flow,
i.e., <code>C</code> and <code>mC_flow</code>.
</li>
<li>
April 28, 2008, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));



 // Medium2.ThermodynamicState state_a2_inflow=
//    Medium2.setState_phX(port_a2.p, inStream(port_a2.h_outflow), inStream(port_a2.Xi_outflow))
//    "state for medium inflowing through port_a2";
//  Medium2.ThermodynamicState state_b2_inflow=
//    Medium2.setState_phX(port_b2.p, inStream(port_b2.h_outflow), inStream(port_b2.Xi_outflow))
 //   "state for medium inflowing through port_b2";//

end PartialFourPortInterfaceThreeMedium;
