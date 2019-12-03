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
AtmosphericCorrection = [(1.01325*ones(1,3)), 0];
Day1Const = [264 183 51 1.437].* ValueCorrection + AtmosphericCorrection;
Day2Const = [268 185 52 1.461].* ValueCorrection + AtmosphericCorrection;
Day3Const = [265 186 51 1.472833].* ValueCorrection + AtmosphericCorrection;
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
state1 = zeros(3,7);
state2 = zeros(3,7);
state2s = zeros(3,7);
state3 = zeros(3,8);
state4 = zeros(3,8);

TCtest = [Averages(:,11) Averages(:,8) Averages(:,9)];
Ptest = [Const(:,3) Const(:,1) Const(:,2)];
stest = [];

for i=1:3
    for j = 1:3
        try
            statest = cell2mat(struct2cell(R22_sh(Ptest(j,i),'T',TCtest(j,i),1)))';
        catch
            statest = 0;
        end

        if statest == 0
            stest{j,i} = 'Saturated';
        else
            stest{j,i} = 'Superheated';
        end
    end
end

for i = 1:3
    state1(i,:) = cell2mat(struct2cell(R22_sh(Const(i,3),'T',Averages(i,11),1)))';
    state2(i,:) = cell2mat(struct2cell(R22_sh(Const(i,1),'T',Averages(i,8),1)))';
    state2s(i,:) = cell2mat(struct2cell(R22_sh(Const(i,1),'s',state1(i,6),1)))';
    state3(i,:) = cell2mat(struct2cell(R22_sat('p', Const(i,2),'x', 0,1)))';
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
bmax = T
beta = (h1-h4)./(h2-h1)
isen_eff_comp = (h2s-h1)./(h2-h1)
power_eff_comp = (mfr.*(h2-h1))/(Voltage*Amps)*1000
%% R22 Heat Calculation
Qin = mfr.*(h1-h4)
Qout = mfr.*(h3-h2)

%% Condenser and Evaporator Efficiencies
QD1 = cell2mat(struct2cell(load('QD1.mat')))';
QD2 = cell2mat(struct2cell(load('QD2.mat')))';
QD3 = cell2mat(struct2cell(load('QD2.mat')))';

QDData = [QD1 ; QD2 ; QD3];
Qd12 = QDData(:,1);
Qd23 = QDData(:,2);

Cond_eff = abs(Qd23./Qout)
Evap_eff = abs((Qin./Qd12).^-1)

%% Statistics
TC1s = Data2(:,5); % [deg C]
RC1s = Data2(:,2); % [%]
Vels = Data2(:,12);% [m/s]
N = length(TC1s);
Ts1 = std(TC1s);
Rs1 = std(RC1s);
Vs1 = std(Vels);
d = 0.5;
invt = tinv(0.95,N-1);
CI_TC1 = tinv(0.95,N-1)*Ts1/sqrt(N);
CI_RC1 = tinv(0.95,N-1)*Rs1/sqrt(N);
CI_Vels = tinv(0.95,N-1)*Vs1/sqrt(N);


N_TC1 = ((invt*Ts1)/d)^2;
N_RC1 = ((invt*Rs1)/d)^2;
N_Vels = ((invt*Vs1)/d)^2;

Nint = 2;
CI = 10000000;
while CI > 0.5002
    nu = Nint;
    Nint = Nint + 1;
    t = tinv(0.95,nu);
    CI = t*Ts1/sqrt(Nint);
end
Nint_TC = Nint;
Nint = 2;
CI = 10000000;
while CI > 0.5002
    nu = Nint;
    Nint = Nint + 1;
    t = tinv(0.95,nu);
    CI = t*Rs1/sqrt(Nint);
end
Nint_RC = Nint;
Nint = 2;
CI = 10000000;
while CI > 0.5002
    nu = Nint;
    Nint = Nint + 1;
    t = tinv(0.95,nu);
    CI = t*Vs1/sqrt(Nint);
end
Nint_Vel = Nint;

TCdiff = [abs(N_TC1 - Nint_TC) 1-abs(N_TC1/Nint_TC)];
RCdiff = [abs(N_RC1 - Nint_RC) 1-abs(N_RC1/Nint_RC)];
Veldiff = [abs(N_Vels - Nint_Vel) abs(N_Vels/Nint_Vel)];

%% Uncertainty
% % Upper Bound Uncertainty
state1 = zeros(3,8);
state2 = zeros(3,7);
state2s = zeros(3,7);
state3 = zeros(3,8);
state4 = zeros(3,8);

Averages2=[Averages(:,1), Averages(:,2), Averages(:,3), Averages(:,4), Averages(:,5)+.5, Averages(:,6)+.5, Averages(:,7)+.5, Averages(:,8)+.5, Averages(:,9)+.5, Averages(:,10)+.5, Averages(:,11)+.5, Averages(:,12), Averages(:,13)];
Const2=[Const(:,1)+0.344738,Const(:,2)+0.344738,Const(:,3)+0.344738,Const(:,4)];

for i = 1:3
    state1(i,:) = cell2mat(struct2cell(R22_sat('T',Averages2(i,11),'x',1,1)))';
    state2(i,:) = cell2mat(struct2cell(R22_sh(Const2(i,1),'T',Averages2(i,8),1)))';
    state2s(i,:) = cell2mat(struct2cell(R22_sh(Const2(i,1),'s',state1(i,6),1)))';
    state3(i,:) = cell2mat(struct2cell(R22_sat('p',Const2(i,2),'x',0,1)))';
    state4(i,:) = cell2mat(struct2cell(R22_sat('T',Averages2(i,10),'h',state3(i,5),1)))';
end

h1up = state1(:,5);
h2up = state2(:,5);
h2sup = state2s(:,5);
h3up = state3(:,5);
h4up = state4(:,5);

v1up = state1(:,3);
v2up = state2(:,3);
v3up = state3(:,3);
v4up = state4(:,3);

% % Lower Bound Uncertainty
state1 = zeros(3,8);
state2 = zeros(3,7);
state2s = zeros(3,7);
state3 = zeros(3,8);
state4 = zeros(3,8);

Averages3=[Averages(:,1), Averages(:,2), Averages(:,3), Averages(:,4), Averages(:,5)-.5, Averages(:,6)-.5, Averages(:,7)-.5, Averages(:,8)-.5, Averages(:,9)-.5, Averages(:,10)-.5, Averages(:,11)-.5, Averages(:,12), Averages(:,13)];
Const3=[Const(:,1)-0.344738,Const(:,2)-0.344738,Const(:,3)-0.344738,Const(:,4)];

for i = 1:3
    state1(i,:) = cell2mat(struct2cell(R22_sat('T',Averages3(i,11),'x',1,1)))';
    state2(i,:) = cell2mat(struct2cell(R22_sh(Const3(i,1),'T',Averages3(i,8),1)))';
    state2s(i,:) = cell2mat(struct2cell(R22_sh(Const3(i,1),'s',state1(i,6),1)))';
    state3(i,:) = cell2mat(struct2cell(R22_sat('p',Const3(i,2),'x',0,1)))';
    state4(i,:) = cell2mat(struct2cell(R22_sat('T',Averages2(i,10),'h',state3(i,5),1)))';
end

h1low = state1(:,5);
h2low = state2(:,5);
h2slow = state2s(:,5);
h3low = state3(:,5);
h4low = state4(:,5);

v1low = state1(:,3);
v2low = state2(:,3);
v3low = state3(:,3);
v4low = state4(:,3);

% % Total Uncertainty
uh1=(h1up-h1low)/2;
uh2=(h2up-h2low)/2;
uh2s=(h2sup-h2slow)/2;
uh3=(h3up-h3low)/2;
uh4=(h4up-h4low)/2;
mfr = Vol_rate./v3;
uV=0.5;
uI=0.05;
uVolrate=0.0005;
uv3=(v3up-v3low)/2;
mfrup=(Vol_rate+uVolrate)./(v3up);
mfrlow=(Vol_rate-uVolrate)./(v3low);
Qinup=mfrup.*(h1up-h4up);
Qinlow=mfrlow.*(h1low-h4low);
Qoutup=mfrup.*(h3up-h2up);
Qoutlow=mfrlow.*(h3low-h2low);

ubeta1a = (h1up-h4)./(h2-h1up)-beta;
ubeta1b = (h1low-h4)./(h2-h1low)-beta;
ubeta2a = (h1-h4up)./(h2-h1)-beta;
ubeta2b = (h1-h4low)./(h2-h1)-beta;
ubeta3a = (h1-h4)./(h2up-h1)-beta;
ubeta3b = (h1-h4)./(h2low-h1)-beta;

ubeta1=(ubeta1a-ubeta1b)/2;
ubeta2=(ubeta2a-ubeta2b)/2;
ubeta3=(ubeta3a-ubeta3b)/2;
ubeta=sqrt(ubeta1.^2+ubeta2.^2+ubeta3.^2);

uisen_eff_comp1a = (h2sup-h1)./(h2-h1)-isen_eff_comp;
uisen_eff_comp1b = (h2slow-h1)./(h2-h1)-isen_eff_comp;
uisen_eff_comp2a = (h2s-h1up)./(h2-h1up)-isen_eff_comp;
uisen_eff_comp2b = (h2s-h1low)./(h2-h1low)-isen_eff_comp;
uisen_eff_comp3a = (h2s-h1)./(h2up-h1)-isen_eff_comp;
uisen_eff_comp3b = (h2s-h1)./(h2low-h1)-isen_eff_comp;

uisen_eff_comp1 = (uisen_eff_comp1a-uisen_eff_comp1b)/2;
uisen_eff_comp2= (uisen_eff_comp2a-uisen_eff_comp2b)/2;
uisen_eff_comp3= (uisen_eff_comp3a-uisen_eff_comp3b)/2;
uisen_eff_comp=sqrt(uisen_eff_comp1.^2+uisen_eff_comp2.^2+uisen_eff_comp3.^2);

power_eff_comp1a = (mfrup.*(h2-h1))/(Voltage*Amps)-power_eff_comp;
power_eff_comp1b = (mfrlow.*(h2-h1))/(Voltage*Amps)-power_eff_comp;
power_eff_comp2a = (mfr.*(h2up-h1))/(Voltage*Amps)-power_eff_comp;
power_eff_comp2b = (mfr.*(h2low-h1))/(Voltage*Amps)-power_eff_comp;
power_eff_comp3a = (mfr.*(h2-h1up))/(Voltage*Amps)-power_eff_comp;
power_eff_comp3b = (mfr.*(h2-h1low))/(Voltage*Amps)-power_eff_comp;
power_eff_comp4a = (mfr.*(h2-h1))/((Voltage+uV)*Amps)-power_eff_comp;
power_eff_comp4b = (mfr.*(h2-h1))/((Voltage-uV)*Amps)-power_eff_comp;
power_eff_comp5a = (mfr.*(h2-h1))/(Voltage*(Amps+uI))-power_eff_comp;
power_eff_comp5b = (mfr.*(h2-h1))/(Voltage*(Amps-uI))-power_eff_comp;

power_eff_comp1 = (power_eff_comp1a-power_eff_comp1b)/2;
power_eff_comp2 = (power_eff_comp2a-power_eff_comp2b)/2;
power_eff_comp3 = (power_eff_comp3a-power_eff_comp3b)/2;
power_eff_comp4 = (power_eff_comp4a-power_eff_comp4b)/2;
power_eff_comp5 = (power_eff_comp5a-power_eff_comp5b)/2;
upower_eff_comp=sqrt(power_eff_comp1.^2+power_eff_comp2.^2+power_eff_comp3.^2+power_eff_comp4.^2+power_eff_comp5.^2);


Cond_eff1a = abs(22.6036./Qout)-Cond_eff; %upper bound Q23 day 2
Cond_eff1b = abs(22.6804./Qout)-Cond_eff; %lower bound Q23 day 2
Cond_eff2a = abs(Qd23./Qoutup)-Cond_eff;
Cond_eff2b = abs(Qd23./Qoutlow)-Cond_eff;

Cond_eff1 = (Cond_eff1a-Cond_eff1b)/2; 
Cond_eff2 = (Cond_eff2a-Cond_eff2b)/2;
uCond_eff = sqrt(Cond_eff1.^2+Cond_eff2.^2);

Evap_eff1a = abs(-13.4341./Qin); %upper bound Q12 day 2
Evap_eff1b = abs(-13.5087./Qin); %lower bound Q12 day 2
Evap_eff2a = abs(Qd12./Qinup);
Evap_eff2b = abs(Qd12./Qinlow);

Evap_eff1 = (Evap_eff1a-Evap_eff1b)/2;
Evap_eff2 = (Evap_eff2a-Evap_eff2b)/2;
uEvap_eff = sqrt(Evap_eff1.^2+Evap_eff2.^2);

uncertainty = {ubeta(2) uisen_eff_comp(2) upower_eff_comp(2) uCond_eff(2) uEvap_eff(2)} %uncertainties using day 2 data
labels = {'B Uncertainty' 'Isen. Compressor N Uncertainty' 'Power Compressor N Uncertainty' 'Condenser N Uncertainty' 'Evaporator N Uncertainty'}
%% Tables
performance = {'' 'B' 'Isentropic Compressor' 'Power Compressor' 'Evaporator' 'Condenser';
                'Day 2' beta(2) isen_eff_comp(2) power_eff_comp(2) Evap_eff(2) Cond_eff(2);
                'Day 3' beta(3) isen_eff_comp(3) power_eff_comp(3) Evap_eff(3) Cond_eff(3)};
            
 stats = {'' 'Number of Measurements';
          'TC1' ceil(N_TC1);
          'RC1' ceil(N_RC1);
          'Velocity 1' ceil(N_Vels)};
      
 uncert = [labels ; uncertainty];
 
%xlswrite('Perf.xlsx',performance,1)
%xlswrite('Stats.xlsx',stats,1)
%xlswrite('Uncert.xlsx',uncert,1)



