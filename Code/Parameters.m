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
a0 = 49;
a1 = 1.94;
a2 = -0.001;
a3 = 3.3e-7;

%% Other parameters
rho = 1.2; %[Kg/m^3]
Beta1 = 0.02;
Beta2 = 0;

%% Flywheel parameters
Jf = 0.026; %[m^2*Kg]

%% GearBox parameteres
tau1 = 1/13.36;
tau2 = 1/10.79;
tau3 = 1/9.15;
tau4 = 1/7.72;
tau5 = 1/6.66;
tau6 = 1/5.76;
tau7 = 1/5;
tau8 = 1/4.35;
tauf = 1/3.73;
Wmax = 1571; %[rad/s]
V1max = tauf*tau1*Rw*Wmax; %[m/s]
V2max = tauf*tau2*Rw*Wmax; %[m/s]
V3max = tauf*tau3*Rw*Wmax; %[m/s]
V4max = tauf*tau4*Rw*Wmax; %[m/s]
V5max = tauf*tau5*Rw*Wmax; %[m/s]
V6max = tauf*tau6*Rw*Wmax; %[m/s]
V7max = tauf*tau7*Rw*Wmax; %[m/s]