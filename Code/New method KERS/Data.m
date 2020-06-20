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
tauk = 1/12;

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
in_cond = 0.00001/(3.6*Rw);

Tmax = 1780; %5*(0.1)^2*Delta(omegaf)/(4*9.55)
Tmax2 = -5504; %600*27.8*0.33 

tau_trasm = 1/(tauk*tauf);