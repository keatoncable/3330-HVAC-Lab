%% HVAC Group Report
clc;clear;close all
%% Import DATA
%Time(1) RH1(2) RH2(3)	RH3(4) Temp1(5) Temp2(6) Temp3(7) TC5(8) TC6(9) TC7(10) TC8(11) Velocity(12) Temperature(13)	
Data1 = load('Data1.txt');
Data2 = load('Data2.txt');
Data3 = load('Data3.txt');
%% Constants
%P1(1) P2(2) P4(3) Flowrate R22(4)   psi    gal/min
ValueCorrection = [(0.00689476 * ones(1,3)*10) (6.306*10^(-5))];
Day1Const = [264 183 51 1.437].* ValueCorrection;
Day2Const = [268 185 52 1.461].* ValueCorrection;
Day3Const = [265 186 51 1.472833].* ValueCorrection;
Const = [Day1Const ; Day2Const ; Day3Const];
Vol_rate = Const(:,4);
Amps = 26.7;
Voltage = 208;
Area = 20 * 23; %inches squared
%% Ideal System Calculations
Averages = zeros(3,length(Data1(1,:)));
for n = 1:length(Data1(1,:))
    Averages(1,n) = mean(Data1(:,n));
    Averages(2,n) = mean(Data2(:,n));
    Averages(3,n) = mean(Data3(:,n));
end

%% Calculating H-S Values
state1 = zeros(3,8);
state2 = zeros(3,7);
state2s = zeros(3,7);
state3 = zeros(3,8);
state4 = zeros(3,8);

for i = 1:3
    state1(i,:) = cell2mat(struct2cell(R22_sat('T',Averages(i,11),'x',1,1)))';
    state2(i,:) = cell2mat(struct2cell(R22_sh(Const(i,1),'T',Averages(i,8),1)))';
    state2s(i,:) = cell2mat(struct2cell(R22_sh(Const(i,1),'s',state1(i,6),1)))';
    state3(i,:) = cell2mat(struct2cell(R22_sat('p',Const(i,2),'x',0,1)))';
    state4(i,:) = cell2mat(struct2cell(R22_sat('T',Averages(i,10),'h',state3(i,5),1)))';
end

h1 = state1(:,5);
h2 = state2(:,5);
h2s = state2s(:,5);
h3 = state3(:,5);
h4 = state4(:,5);

v1 = state1(:,3);
v2 = state2(:,3);
v3 = state3(:,3);
v4 = state4(:,3);

%% Efficiencies 
mfr = Vol_rate./v3;
beta = (h1-h4)./(h2-h1);
isen_eff_comp = (h2s-h1)./(h2-h1);
power_eff_comp = (mfr.*(h2-h1))/(Voltage*Amps);

%% R22 Heat Calculation
Qin = mfr.*(h1-h4);
Qout = mfr.*(h2-h3);

test = 'test';
test2 = 'hey';
%% Statisitcs

