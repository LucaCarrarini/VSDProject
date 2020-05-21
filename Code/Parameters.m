%% Vehicle parameters
M = 700; %[kg]
Cx = 0.9;
Sf = 1.5; %[m^2]

%% Wheels parameters
Rw = 0.33; %[m]
Mw = 9; %[kg]
Jw = 1/2*(Rw^2*Mw); %[m^2*Kg]

%% Engine parameters
Je = 0.3; %[m^2*Kg]

%% Other parameters
rho = 1.2; %[Kg/m^3]
Beta1 = 0.02;
Beta2 = 0;

%% Flywheel parameters
Jf = 0.026; %[m^2*Kg]

%% GearBox parameters

%These data are taken from La Ferrari datasheet
%7 ratios transmission

tau1 = 1/3.08;
tau2 = 1/2.18;
tau3 = 1/1.63;
tau4 = 1/1.29;
tau5 = 1/1.03;
tau6 = 1/0.84;
tau7 = 1/0.69;
tauf = 1/4.38;
Wmax = 942.4; %[rad/s]
V1max = tauf*tau1*Rw*Wmax; %[m/s]
V2max = tauf*tau2*Rw*Wmax; %[m/s]
V3max = tauf*tau3*Rw*Wmax; %[m/s]
V4max = tauf*tau4*Rw*Wmax; %[m/s]
V5max = tauf*tau5*Rw*Wmax; %[m/s]
V6max = tauf*tau6*Rw*Wmax; %[m/s]
V7max = tauf*tau7*Rw*Wmax; %[m/s]