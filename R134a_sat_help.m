function R134a = R134a_sat(name1, value1, name2, value2, units) 
% -------------------------------------------------------------------------
% R134a_sat.m
% -------------------------------------------------------------------------
% Fundamentals of Thermodynamics
% A function that uses two independent properties of saturated R134a 
% (as within the vapor dome) and defines as the output a structure with
% all relevant properties at the same state. Since the state is within the
% vapor dome, a combination of temperature or pressure, along with any 
% other valid property will fix the state. 
% -------------------------------------------------------------------------
% Syntax: R134a = R134a_sat(name1, value1, name2, value2, units) 
% Inputs: name1 = the name of the first property ('T' or 'p')
%         value1 = the value of the first property
%         name2 = the name of the second property ('v','u','h','s' or 'x')
%         value2 = the value of the second property
%         units = 1 (SI units) or 2 (English units)
% Output: R134a = the structure with all relevant properties in relevant
%         fields (T, p, v, u, h, s, x and message).  
%         message is 1 if state is in vapor dome, 0 if not.
%
% Valid properties and units:
%         T - temperature (C or F)
%         p - pressure (bar or lbf/in^2)
%         v - specific volume (m^3/kg or ft^3/lb)
%         h - specific enthalpy (kJ/kg or Btu/lb)
%         u - specific internal energy (kJ/kg or Btu/lb)
%         s - specific entropy (kJ/kg*K or Btu/lb*R)
%         x - quality (mass of vapor as a unitless fraction between 0-1)
% -------------------------------------------------------------------------
% Notes: Use consistent SI or English units.
% All specific properties defined per unit mass. 
%--------------------------------------------------------------------------
% Example
% R134a = R134a_sat('T', 0, 'v', 0.01, 1)
% This example corresponds to saturated refrigerant 134a at 0 degrees C, 
% and with a specific volume of 0.01 m^3/kg
%--------------------------------------------------------------------------
% Simon Ruiz, Undergraduate student, Mechanical Engineering 
% Dr.Priya Goeser, Professor
% Original: March 11th, 2013
% Modified: July 21, 2016
% Copyright (c) 2013 Simon Ruiz, All Rights Reserved.
% Engineering Studies Program, Armstrong State University. 
% -------------------------------------------------------------------------