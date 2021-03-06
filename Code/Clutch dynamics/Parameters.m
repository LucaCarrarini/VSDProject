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
w_fw_max = 6283; %[rad/s]
rho_fw = 2267; %[kg/m^3]
h = 0.1; %[m]
r_o = 0.1; %[m]
r_i = 0.05; %[m]
% SOC_max = 0.25*pi*h*rho_mat*(r_o^4-r_i^4)

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


%% Dynamic clutch engagement

Jeq = 0.5*(Je/(tauf)+M*Rw^2+4*Jw);
Jc = 0.01;
b_cv = 0.03; b_ck = 0.03;
mi_d = 0.4;
k = 0.098;
FnMAX = 5000;

%tau1 = 1/13.36;
%tau2 = 1/10.79;
%tau3 = 1/9.15;
%tau4 = 1/7.72;
%tau5 = 1/6.66;
%tau6 = 1/5.76;
%tau7 = 1/5;
%tau8 = 1/4.35;
%tauf = 1/3.73;
Wmax = 942.4; %[rad/s]
V1max = tauf*tau1*Rw*Wmax; %[m/s]
V2max = tauf*tau2*Rw*Wmax; %[m/s]
V3max = tauf*tau3*Rw*Wmax; %[m/s]
V4max = tauf*tau4*Rw*Wmax; %[m/s]
V5max = tauf*tau5*Rw*Wmax; %[m/s]
V6max = tauf*tau6*Rw*Wmax; %[m/s]
V7max = tauf*tau7*Rw*Wmax; %[m/s]


%% LQR control

%state: omega_cv, omega_cv-omega_ck, Fn
A = [-b_cv/Jeq, 0, -k/Jeq;...
    -b_cv/Jeq+b_ck/Jc, -b_ck/Jc, -(k/Jeq+k/Jc);...
    0 0 0];

B = [0;0;1];

% Initial condition
x_0 = [210;210;0]; % [rad/s, rad/s, N]

Q =[0, 0,    0;...
    0,  10,   0;...
    0,  0,    0];

R = 1;

%[K,P,e] = lqr(A,B,Q,R);

% Disturbance vector field
% Bd =  [1/Jc, 0;...
%        1/Jc, 1/Jeq;...
%        0,    0];
%

Bd =  [1/Jeq;
       1/Jeq; 
       0];    
   
% Disturbance constant
T_eq = -3970; % [Nm]
%T_l = 0; % [Nm]
d = [T_eq];

Gamma = Bd*d;

%% Compute control
control.P = flip(out.P,3);
control.M = flip(out.M,1);
control.h = flip(out.h,1);
control.U = flip(out.U,3);
control.F = flip(out.F,1);
control.K = flip(out.K,1);
v = (control.F(1))\(-control.U(:,:,1)*x_0 - control.K(1));

sigP.time = out.tout;
sigP.signals.values = control.P;
sigP.signals.dimensions = [3,3];

sigM.time = out.tout;
sigM.signals.values = (control.M);
sigM.signals.dimensions = 3;

sigh.time = out.tout;
sigh.signals.values = (control.h);
sigh.signals.dimensions = 3;