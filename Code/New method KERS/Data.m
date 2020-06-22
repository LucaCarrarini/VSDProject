%% Parametri progettuali
clc
clear all
close all
%% Wheels
Rw = 0.33; 
Mw = 9; 
Jw = 1/2*(Rw^2*Mw);
%% Engine
Je = 0.3;
%% vehicle
M = 700+70; % Vehicle+Pilot 
Jeq_v = M*Rw^2; % neglecting engine and wheels inertias
%% Flywheel
Rf = 0.1;
Mf = 5;
Jf = 1/2*(Rf^2*Mf);
wf_max = 64500/9.55;
%% Aerodynamic drag
Cx = 0.9;
Sf = 1.5;
rho = 1.2;
d = 0.5*rho*Cx*Sf;
%% Torotrak
tau_min = 0.4; 
tau_max = 2.5;
%% Gear Box
tau1 = 1/3.08;
tau2 = 1/2.18;
tau3 = 1/1.63;
tau4 = 1/1.29;
tau5 = 1/1.03;
tau6 = 1/0.84;
tau7 = 1/0.69;
tauf = 1/4.38;
tauk_charge = 1/12;
tauk_discharge = 11*tauk_charge;

Wmax = 942.4; %[rad/s]
V1max = tauf*tau1*Rw*Wmax; %[m/s]
V2max = tauf*tau2*Rw*Wmax; %[m/s]
V3max = tauf*tau3*Rw*Wmax; %[m/s]
V4max = tauf*tau4*Rw*Wmax; %[m/s]
V5max = tauf*tau5*Rw*Wmax; %[m/s]
V6max = tauf*tau6*Rw*Wmax; %[m/s]
V7max = tauf*tau7*Rw*Wmax; %[m/s]

% Efficiency
eta = 1;

%% KERS clutch

R_clutch = 0.18;
mi_d = 0.4;
k = 4/3*mi_d*R_clutch;
Fn_max = 10000; % [N]

% Time at which the kers discharge starts
time_push = 25; % [sec]

Tmax = 240; %5*(0.1)^2*Delta(omegaf)/(4*9.55)
Tmax2 = -2457; %700*9.8*0.33 

tau_trasm_charge = 1/(tauk_charge*tauf);
%tau_trasm_discharge = tau_trasm_charge/10;
tau_trasm_discharge = 1/(tauk_discharge*tauf);