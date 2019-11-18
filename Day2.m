clear
clc
close all

D = load('Data2.txt');
%%
A = .0254*20*.0254*23; %inlet area [m^2]
Vel = mean(D(:,12)); %avg inlet velocity [m/s]
fr = A*Vel; %[m^3/s]
T1 = mean(D(:,5)) + 273; %[K]
H1 = mean(D(:,2)); % [%]
patm = 101325; % [Pa]
Ra = 286.9; % [J/kg K]
Rw = 461.5; % [J/kg K]
pga = 2339;
pgb = 2487;
Ta = floor(T1);
Tb = ceil(T1);
pg = (pgb-pga)/(Tb-Ta)*(T1-Ta)+pga; % saturation pressure
pv = H1*pg/100; % partial vapor pressure
w1 = 0.622*pv/(patm-pv); % humidity ratio
rho = patm/(Ra*T1)*(1-0.378*pv/patm); % [kg/m^3]
mdot = rho*fr; % mass flow rate [kg/s]
mdota = mdot/(w1+1); % mass flow rate of dry air which stays constant throughout process
mdotv1 = mdot-mdota; % mass flow rate of water vapor which is changed by dehumidification process

%%
T2 = mean(D(:,6)) + 273;
H2 = mean(D(:,3));
T3 = mean(D(:,7)) + 273;
H3 = mean(D(:,4));
pga = 3169;
pgb = 3363;
Ta = floor(T3);
Tb = ceil(T3);
pg = (pgb-pga)/(Tb-Ta)*(T3-Ta)+pga; % saturation pressure
pv = H3*pg/100; % partial vapor pressure
% pnew = T3*patm/T1;
w3 = 0.622*pv/(patm-pv); % humidity ratio
w2 = w3; %this is a known/assumed relationship
mdotw = mdota*(w1-w2);
mdotv2 = mdota*w2;

%%
ha = 290.16;
hb = 295.17;
Ta = 290;
Tb = 295;
ha1 = (hb-ha)/(Tb-Ta)*(T1-Ta)+ha;
ha = 2538.1;
hb = 2539.9;
Ta = 293;
Tb = 294;
hv1 = (hb-ha)/(Tb-Ta)*(T1-Ta)+ha;
ha = 285.14;
hb = 290.16;
Ta = 285;
Tb = 290;
ha2 = (hb-ha)/(Tb-Ta)*(T2-Ta)+ha;
ha = 2525.3;
hb = 2527.1;
Ta = 286;
Tb = 287;
hv2 = (hb-ha)/(Tb-Ta)*(T2-Ta)+ha;
ha = 54.6;
hb = 58.8;
hw = (hb-ha)/(Tb-Ta)*(T2-Ta)+ha;
ha = 295.17;
hb = 300.19;
Ta = 295;
Tb = 300;
ha3 = (hb-ha)/(Tb-Ta)*(T3-Ta)+ha;
ha = 2547.2;
hb = 2549;
Ta = 298;
Tb = 299;
hv3 = (hb-ha)/(Tb-Ta)*(T3-Ta)+ha;

Qdot12 = mdota*((ha2+w2*hv2)-(ha1+w1*hv1)+(w1-w2)*hw);
Qdot23 = mdota*(ha3-ha2) + mdotv2*(hv3-hv2);