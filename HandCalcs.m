clear
clc
close all

set(0,'defaulttextInterpreter','latex') %latex axis labels
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
format shortEng

%--------------------------------------------------------------------------
% LOADS
%--------------------------------------------------------------------------

%variables

air_compressor_omega = 55 * (2 * pi) / 60;
generator_omega = 450 * (2 * pi) / 60;
gen_diameter = 500e-3;
air_compressor_power = 100e+3;
generator_power = 2e+3;
beta = 12; %degrees
Sy = 410e6; % [Pa]
UTS = 680e6; % [Pa]
SCF = 1.6; % Factor of safety for the shaft
friction_factor = 0.3; % Steel on cast iron
shaftOuter_diameter = 0.110; % Shaft outer diameter [m]
shaftInner_diameter = 0.100; % Shaft inner diameter [m]
clutchOuter_diameter = 0.45; % Clutch inner diameter [m]
%distances
distance_between_bearings = 900e-3;
distance_COM_2_outer_bearing = 650e-3;
distance_belt_2_inner_bearing = 320e-3;
distance_driveshaft_2_inner_brearing = distance_belt_2_inner_bearing + 50e-3;
PCD= 500e-3;
weight_force = 850 * 9.81;
axial_force = 18e3; % force [N]


%belt force calcs
T_belt = torqueShaft(generator_power,generator_omega); % Torque from gen
belt_force = beltForce(T_belt,PCD/2);

% bearing force
T_comp = torqueShaft(air_compressor_power,air_compressor_omega);
[A_b,B_b] = bearings(T_comp, beta,distance_between_bearings);

B = [1 1;...
    distance_COM_2_outer_bearing distance_between_bearings + distance_COM_2_outer_bearing];

Reactions = linsolve(B,[-belt_force + weight_force;...
    -belt_force * (distance_driveshaft_2_inner_brearing + distance_between_bearings + distance_COM_2_outer_bearing)]); 
% moment calcs
A_max = 14420; %A_b + Reactions(1);
B_max = -6251; %B_b + Reactions(2);

T_total = T_comp - T_belt; % Torque of belt takes away from torque of main shaft
%T_com = torqueShaft(airCompressorPower,airCompressorOmega)

xArray = linspace(0, (650 + 900 + 320 + 50) * 10 ^ -3, 1000);
[vArray,mArray] = momentDiagram(xArray, ...
    -850 * 9.81, 0, ...
    A_max,650e-3, ...
    B_max,(650+900) * 10 ^ -3, ...
    belt_force, (650 + 900 + 320) * 10 ^ -3);
plot(xArray,vArray,'lineWidth',3)
ylabel("$V$","fontSize",14)
xlabel("$x$","fontSize",14)
figure
plot(xArray,mArray,'lineWidth',3)
ylabel("$M$","fontSize",14)
xlabel("$x$","fontSize",14)

max_moment = max(-1 * mArray);

min_shaft_diameter = minShaftDiameter(T_comp, max_moment);

FOS = ShaftFOS(T_total, max_moment, Sy, SCF, shaftOuter_diameter);

key_width = 0.25 * shaftOuter_diameter; % Width of key (Industry standard to choose 1/4 diamter of shaft 
key_depth = 0.125 * shaftOuter_diameter; % Depth of shaft (1/2 way above shaft)

% Display Results

key_length_shear = keyLengthShear(T_total, key_width, Sy, shaftOuter_diameter);

key_length_compressive = keyLengthCompressive(T_total, key_depth, Sy, shaftOuter_diameter);

% Shoulder calculations

shoulder_stress = shoulderCompressiveStress(shaftOuter_diameter, shaftInner_diameter, axial_force);

% Clutch calculations

actuating_force = UniformWear(T_comp, friction_factor, clutchOuter_diameter, shaftInner_diameter)



% Belt
fprintf("Belt force [N]: %2.2f, Tension of Belt: %2.2f\n", belt_force, T_belt);

% Shaft Diamter
fprintf("Minimum shaft diameter: %2.4f mm\n",min_shaft_diameter*1e3);
fprintf("Max Moment: %2.4f Nm\n",max_moment);
fprintf("Shaft Diameter with FOS(%2.2f): %2.2f mm \n", FOS, shaftOuter_diameter*1000);

% Keyway results
fprintf("Length of key under shear stress: %2.2f mm\n", key_length_shear*1000);
fprintf("Length of key under compressive stress: %2.2f mm \n", key_length_compressive*1000);


% Shoulder results
fprintf("Compressive stress on shoulder is: %2.2f MPa \n", shoulder_stress/1e6);

%--------------------------------------------------------------------------
% IMPORTANT EQUATIONS
%--------------------------------------------------------------------------
function [shoulder_comp_stress] = shoulderCompressiveStress(shaftOuter_diameter, shaftInner_diameter, axial_force)
    
    surface_area = (pi/4) * (shaftOuter_diameter^2 - shaftInner_diameter^2);
    shoulder_comp_stress = axial_force/surface_area;
end

function [key_length_comp] = keyLengthCompressive(T_total, key_depth, Sy, shaft_diameter)
    % Calculate the length of the key required to not fail under
    % compressive loads
    
    key_length_comp = (4*T_total)/(shaft_diameter*Sy*key_depth);
    
end

function [key_length_shear] = keyLengthShear(T_total, key_width, Sy, shaft_diameter)
    % Calculate the length of the 
    key_length_shear = (2*T_total)/((Sy/2) * key_width * shaft_diameter);
end

function[FOS] = ShaftFOS(T_total, max_moment, Sy, SCF, shaft_diameter)
    % Calculate the shaft diameter with accounting for a specific FOS
    max_shear = 16/(pi*shaft_diameter^3) * sqrt(max_moment^2 + T_total^2);
    tensile_stress = 2 * max_shear;
    FOS = Sy/(tensile_stress*SCF);
    
end

function[d] = minShaftDiameter(torque, bending_moment)
    % Constants
    Cm = 1.5;
    Ct = 1.0;
    max_shear = 0.18 * 680e6; % UTS*0.18
    
    
    d = ( (5.1/max_shear) * ( (Cm*bending_moment)^2 + (Ct*torque)^2)^0.5)^(1/3);

end


function[V,M] = momentDiagram(x,f1,d1,f2,d2,f3,d3,f4,d4)
    V = heaviside(x - d1) * f1 ...
        + heaviside(x - d2) * f2 ...
        + heaviside(x - d3) * f3 ...
        + heaviside(x - d4) * f4;
    M = heaviside(x - d1) .* f1 .* (x)...
        + heaviside(x - d2) .* f2 .* (x-d2)...
        + heaviside(x - d3) .* f3 .* (x-d3)...
        + heaviside(x - d4) .* f4 .* (x-d4);
    
end

function [T] = torqueShaft(P,omega)
    T = P / omega;
end

function [F] = UniformWear(T_W,f,D,d)
    F = (4*T_W)/(f*(D+d));
end

function [T_W] = UniformPressure(F,f,D,d)
    T_W = (F * f) / (3) * ( (D^3 - d ^3) / (D^2 - d^2) );
end

function [A,B] = bearings(T_DI, beta,a)
    A = -(T_DI * tand(beta)) / (a); %degrees
    B = +(T_DI * tand(beta)) / (a); %degrees
end

function [F] = beltForce(T,R)
    F = (T) / (R); % Assumed  belt is like chain (one side is zero tension like chain)
end
