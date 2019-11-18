clear
clc
close all

D = load('D1.txt');
%%
A = .0254*20*.0254*23; %inlet area [m^2]
Vel = mean(D(:,12)); %avg inlet velocity [m/s]
fr = A*Vel; %[m^3/s]
T1 = mean(D(:,5)) + 273; %[K]
H1 = mean(D(:,2)); % [%]
patm = 101325; % [Pa]
Ra = 286.9; % [J/kg K]
Rw = 461.5; % [J/kg K]
pga = 2645;
pgb = 2810;
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
pga = 3363;
pgb = 3567;
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
ha = 295.17;
hb = 300.19;
Ta = 295;
Tb = 300;
ha1 = (hb-ha)/(Tb-Ta)*(T1-Ta)+ha;
ha3 = (hb-ha)/(Tb-Ta)*(T3-Ta)+ha;
ha = 2541.7;
hb = 2543.5;
Ta = 295;
Tb = 296;
hv1 = (hb-ha)/(Tb-Ta)*(T1-Ta)+ha;
ha = 280.13;
hb = 285.14;
Ta = 280;
Tb = 285;
ha2 = (hb-ha)/(Tb-Ta)*(T2-Ta)+ha;
ha = 2521.6;
hb = 2523.4;
Ta = 284;
Tb = 285;
hv2 = (hb-ha)/(Tb-Ta)*(T2-Ta)+ha;
ha = 46.2;
hb = 50.41;
hw = (hb-ha)/(Tb-Ta)*(T2-Ta)+ha;
ha = 2549;
hb = 2550.8;
Ta = 299;
Tb = 300;
hv3 = (hb-ha)/(Tb-Ta)*(T3-Ta)+ha;

Qdot12 = mdota*((ha2+w2*hv2)-(ha1+w1*hv1)+(w1-w2)*hw);
Qdot23 = mdota*(ha3-ha2) + mdotv2*(hv3-hv2);