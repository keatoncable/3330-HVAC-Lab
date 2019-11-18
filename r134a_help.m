function R134a = R134a_sh(pressure, name2, value2, units) 
% -------------------------------------------------------------------------
% R134a_sh.m
% -------------------------------------------------------------------------
% Fundamentals of Thermodynamics
% A function that uses two independent properties of superheated 
% refrigerant 134a vapor (steam) and defines as the output a structure with
% all relevant properties at the same state. Since the state is in the 
% superheated region, a combination of pressure, along with any 
% other valid property will fix the state.  
% -------------------------------------------------------------------------
% Syntax: R134a = R134a_sh(pressure, name2, value2, units) 
% Inputs: pressure= the value of pressure at given state
%         name2 = the name of the second property ('T','v','u','h' or 's')
%         value2 = the value of the second property
%         units = 1 (SI units) or 2 (English units)
% Output: R134a = the structure with all relevant properties in relevant
%         fields (T, p, v, u, h, s and x).  Though undefined (NaN), x is
%         also included to be consistent with the saturated functions. 
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
% R134a = R134a_sh(10, 'T', 120, 1)
% This example corresponds to superheated R134a at a pressure of 
% 10 bar, and with a temperature of 120 degrees C.
%--------------------------------------------------------------------------
% Simon Ruiz, Undergraduate student, Mechanical Engineering 
% Dr.Priya Goeser, Professor
% Original: March 12th, 2013
% Modified: July 21, 2016
% Copyright (c) 2013 Simon Ruiz, All Rights Reserved.
% Engineering Studies Program, Armstrong State University. 
% -------------------------------------------------------------------------