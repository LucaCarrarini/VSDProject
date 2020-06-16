%% Parametri progettuali

clc
clear all
close all

m = 1343; 
R = 0.31;
Jeq_v = 120;
Jf = 0.026;
%0.5rho*Cx*Sf
d = 0.49;
wf_max = 60000/9.55;

% Torotrak
tau_min = 0.4; tau_max = 2.5;

% Efficiency
eta = 1;

Tmax = 240;
Tmax2 = -3675;

tau_trasm = 62.8;
tau_trasm_2 = 0.194;