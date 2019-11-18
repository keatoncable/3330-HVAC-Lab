function [V] = R22State(Value1,Value2,Parameter1,Parameter2,State)
%Temp [deg C] (1)
%Pressure [MPa] (2)
%Specific Volume [m^3/kg] (3)
%Enthalpy [kJ/kg] (4)
%Entropy [kJ/kg K] (5)
%Location is the vector of which values are val1 and val2
%State includes 'L','V', or quality
load R22.mat;
Table = R22_Values;
Table(:,3) = 1./R22_Values(:,3);
V = zeros(1,6);
% Find Location and State
Parameters = ['t','p','v','h','s'];
for n = 1:5
    if strcmpi(Parameter1,Parameters(n))
        Location1 = n;
    end
    
    if strcmpi(Parameter2,Parameters(n))
        Location2 = n;
    end
end

Min = Table(:,Location1) - Value1;
idx = find(Min <= 0, 1,'last');


    
if Location1 + Location2 == 3
    
    if strcmpi(State,'f')
        Quality = 0;
    elseif strcmpi(State,'v')
        Quality = 1;
    elseif State <= 1
        Quality = State;
    end
    
    V(Location1) = Value1;
    V(Location2) = Value2;
    V(6) = Quality;
    
    vf = (Table(idx+1,3)-Table(idx,3))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,3);
    vv = (Table(idx+1,4)-Table(idx,4))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,4);
    V(3) = vf + Quality*(vv-vf);
    
    hf = (Table(idx+1,5)-Table(idx,5))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,5);
    hv = (Table(idx+1,6)-Table(idx,6))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,6);
    V(4) = hf + Quality*(hv-hf);
    
    sf = (Table(idx+1,7)-Table(idx,7))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,7);
    sv = (Table(idx+1,8)-Table(idx,8))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,8);
    V(5) = sf + Quality*(sv-sf);
else
    if Location2 <= 2
        save1 = Value1;
        save2 = Location1;
        Value1 = Value2;
        Location1 = Location2;
        Value2 = save1;
        Location2 = save2;
    end

    if Location2 == 4
        Hf = (Table(idx+1,5)-Table(idx,5))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,5);
        Hv = (Table(idx+1,6)-Table(idx,6))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,6);
    elseif Location2 == 5
        Sf = (Table(idx+1,7)-Table(idx,7))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,7);
        Sv = (Table(idx+1,8)-Table(idx,8))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,8);
    end

    if strcmpi(State,'f')
        Quality = 0;
    elseif strcmpi(State,'v')
        Quality = 1;
    elseif State <= 1
        Quality = State;
    else
        if Location2 == 4
            Quality = (Value2 - Hf)/(Hv-Hf);
        elseif Location2 == 5
            Quality = (Value2 - Sf)/(Sv-Sf);
        end
    end
    
    V(Location1) = Value1;
    V(Location2) = Value2;
    V(6) = Quality;

    for n = 1:4
        if n <=2
            if V(n) == 0
                V(n) = (Table(idx+1,n)-Table(idx,n))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,n);
            end
        elseif n > 2
            if V(n) ==0 && n ==3
                vf = (Table(idx+1,3)-Table(idx,3))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,3);
                vv = (Table(idx+1,4)-Table(idx,4))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,4);
                V(n) = vf + Quality*(vv-vf);
            elseif V(n) == 0 && n == 4
                hf = (Table(idx+1,5)-Table(idx,5))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,5);
                hv = (Table(idx+1,6)-Table(idx,6))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,6);
                V(n) = hf + Quality*(hv-hf);
            elseif V(n) == 0 && n ==5
                sf = (Table(idx+1,7)-Table(idx,7))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,7);
                sv = (Table(idx+1,8)-Table(idx,8))/(Table(idx+1,Location1)-Table(idx,Location1))*(Value1-Table(idx,Location1)) + Table(idx,8);
                V(n) = sf + Quality*(sv-sf);
            end
        end
    end
end


end