clear
clc
close all

A1 = .0254*20*.0254*23; %inlet area [m^2]

T1 = (80-32)*5/9 + 273;
T3 = T1;
H1 = .65;
H2 = .8;
H3 = .4;
patm = 101325;

Ra = 286.9; % [J/kg K]
Rw = 461.5; % [J/kg K]
pga = 3363;
pgb = 3567;
pg13 = (pgb-pga)/(300-299)*(T1-299)+pga; % saturation pressure
pv1 = H1*pg13;
pv3 = H3*pg13;
w1 = 0.622*pv1/(patm-pv1); %humidity ratio in
w3 = 0.622*pv3/(patm-pv3); %humidity ratio out
w2 = w3;
pv2 = pv3;
pg2 = pv2/H2;
T2 = (pg2-1705)/(1818-1705)+15 + 273;
rho1 = patm/(Ra*T1)*(1-0.378*pv1/patm); % [kg/m^3]

ha1 = (300.19-295.17)*(T1-295)/5 + 295.17;
ha2 = (290.16-285.14)*(T2-285)/5 + 285.14;
hv1 = (2550.8-2549)*(T1-299)+2549;
hv2 = (2530.8-2528.9)*(T2 - 273-15)+2528.9;
hw = (67.19-62.99)*(T2 - 273-15)+62.99;

Qdot12 = -15.1777;
mdota = (ha2+w2*hv2-ha1-w1*hv1+(w1-w2)*hw)/Qdot12;
mdot = mdota*(w1+1);
v = mdot/(rho1*A1);