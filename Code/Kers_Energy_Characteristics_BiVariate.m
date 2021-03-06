clc
clear all
close all
%% Vehicle parameters
M = 700; %[kg]
%% Wheels parameters
Rw = 0.33; %[m]
Mw = 9; %[kg]
Jw = 1/2*(Rw^2*Mw); %[m^2*Kg]
%% Engine parameters
Je = 0.3; %[m^2*Kg]
Wmax = 942.4; %[rad/s]
%% Flywheel parameters
Jf = 0.026; %[m^2*Kg]
taukmin = 1/12;
%% Clutch parameters
Jc = 0.01; %[m^2*Kg]
%% CVT parameters
taumin = 0.4; 
taumax = 2.5;
%% GearBox parameters
tau1 = 1/3.08;
tau2 = 1/2.18;
tau3 = 1/1.63;
tau4 = 1/1.29;
tau5 = 1/1.03;
tau6 = 1/0.84;
tau7 = 1/0.69;
tauf = 1/4.38;
%% Gear shift strategy
V1max = tauf*tau1*Rw*Wmax; %[m/s]
V2max = tauf*tau2*Rw*Wmax; %[m/s]
V3max = tauf*tau3*Rw*Wmax; %[m/s]
V4max = tauf*tau4*Rw*Wmax; %[m/s]
V5max = tauf*tau5*Rw*Wmax; %[m/s]
V6max = tauf*tau6*Rw*Wmax; %[m/s]
% %% Energy transfer acceleration (Bi-Variate)
% 
% %Generate array with flywheel test speeds
% %omegaf = linspace(64500*0.3/9.55, 64500/9.55, 100);
% 
% % Initial velocity definition
% V0 = (0/3.6);
% omegaf = 64.5*10^3/9.55;
% %Generate array with vehicle test speeds
% % Initial velocity = 0 m/s
% %v = zeros(1,100); %[m/s]
% 
% %Generate array with test ratios
% tauCVT = linspace(taumin, taumax, 10);
% tauk = linspace(taukmin, taukmin*11, 8);
% tauCVT0 = taumin;
% 
% % for i=1:8
% %     tauk(i) = (log(i+exp(1)))*taukmin;
% % end
% 
% if V0<V1max
%    taugb = tau1;
% end    
% if V0>V1max && V0<V2max
%    taugb = tau2;
% end    
% if V0>V2max && V0<V3max
%    taugb = tau3;
% end    
% if V0>V3max && V0<V4max
%    taugb = tau4;
% end    
% if V0>V4max && V0<V5max
%    taugb = tau5;
% end    
% if V0>V5max && V0<V6max
%    taugb = tau6;
% end    
% if V0>V6max
%    taugb = tau7;   
% end
% %Generate energy curves
% E = ones(length(tauk), length(tauCVT));
% %E0 = ones(1, length(v));
% for i=1:length(tauk)
%     for j=1:length(tauCVT)
%         taugb = 1;
%         %Jeqv = M*((tauf*taugb*Rw/tauCVT(i))^2)+Je/(tauCVT(i))^2+4*Jw*(tauf*taugb/tauCVT(i))^2+Jc;
%         Jeqv = M*((tauf*taugb*Rw*tauCVT0)^2)+Je*(tauCVT0)^2+4*Jw*(tauf*taugb*tauCVT0)^2+Jc;
%         Jeqk = Jf/(taukmin)^2+Jc;
%         %E0kers = 0.5*Jf*omegaf(j)^2;
%         E0kers = 0.5*Jf*omegaf^2;
%         DeltaE = E0kers*M/Jf*(Jeqk/(Jeqv+Jeqk))^2*(tauf*taugb*Rw*tauCVT(j)*tauk(i))^2;
%         E(i,j) = DeltaE;
%         %E0(j) = E0kers;
%     end
% end
% 
% figure()
% for i=1:length(tauk)
%     plot(tauCVT, E(i,:));
%     %text(E0(100), E(i, 100), num2str(gamma(i)));
%     hold on;
% end
% %omega_f_min_rpm = 64500*0.3/9.55; %[rad(sec]
% %Energy_min = 1/2*Jf*(omega_f_min_rpm)^2;
% %xline(Energy_min, 'black--');
% title('Energy transfer during acceleration');
% xlabel('\tau_{CVT}');
% ylabel('\DeltaE_{Vehicle} [J]');
% xlim([taumin, taumax]);
% %ylim([-300 4000])
% legend({num2str(tauk')});
% 
% figure()
% colormap('jet')
% [X, Y] = meshgrid(tauCVT, tauk);
% %psi = Jeq/(Jf)*Jc^2/(Jc+Jeq)^2*((Y./X).^2);
% %rho = (Jeq^2*Jc)/((Jc+Jeq)^2)*Y./X;
% Z = E0kers*M/Jf*(Jeqk/(Jeqv+Jeqk))^2*(tauf*taugb*Rw*X.*Y).^2;
% surf(X,Y,Z)
% xlabel('\tau_{CVT}')
% ylabel('\tau_k')
% title('Energy transfer with variables \tau_k,\tau_{CVT}')
% zlabel('\DeltaE_{Vehicle} [J]')
% 
% figure()
% for i=1:length(tauk)
%     DeltaV = sqrt(2*E(i,:).*sign(E(i,:))*M).*sign(E(i,:));
%     plot(tauCVT, DeltaV*3.6);
% %   text(omegaf(100)*9.55, vf(100)*3.6, num2str(gamma(i)));
%     hold on;
% end
% 
% %omega_f_min_rpm = 64500*0.3; %[rpm]
% %xline(omega_f_min_rpm, 'black--');
% title('Velocity comparison during acceleration');
% xlabel('\tau_{CVT}');
% ylabel('\DeltaV_{f,Vehicle} [Km/h]');
% xlim([taumin, taumax]);
% %ylim([-5 30]);
% legend({num2str(tauk')});

%% Energy transfer deceleration (bi-variate)

V0 = 250/3.6;
tauCVT = linspace(taumin, taumax, 100);
tauk = linspace(taukmin, taukmin*11, 8);
tauCVT0 = taumax;

% for i=1:8
%     tauk(i) = (log(i+exp(1)))*taukmin;
% end


%Generate energy curves
E = ones(length(tauk),length(tauCVT));
%E0 = ones(1, length(tauCVT));
if V0<V1max
   taugb = tau1;
end    
if V0>V1max && V0<V2max
   taugb = tau2;
end    
if V0>V2max && V0<V3max
   taugb = tau3;
end    
if V0>V3max && V0<V4max
   taugb = tau4;
end    
if V0>V4max && V0<V5max
   taugb = tau5;
end    
if V0>V5max && V0<V6max
   taugb = tau6;
end    
if V0>V6max
   taugb = tau7;   
end
for i=1:length(tauk)
    for j=1:length(tauCVT)
        taugb = 1;
        v(j) = 250/3.6;
        Jeqv = M*((tauf*taugb*Rw*tauCVT0)^2)+Je*(tauCVT0)^2+4*Jw*(tauf*taugb*tauCVT0)^2+Jc;
        Jeqk = Jf/(taukmin*11)^2+Jc;
        T0 = 0.5*M*v(j)^2;
        DeltaT = T0*Jf/M*(Jeqv/(Jeqv+Jeqk))^2*(1/(tauf*taugb*Rw*tauk(i)*tauCVT(j)))^2;
        E(i,j) = DeltaT;
        E0(j) = T0;
    end
end

figure()
for i=1:length(tauk)
    plot(tauCVT, E(i,:));
    %text(tauCVT(100), E(i, 100), num2str(tauk(i)));
    hold on;
end
plot(tauCVT, 593*10^3*ones(1, 100), 'black--');
ylim([0, 2.5*10^5]);
xlim([taumin, taumax]);
title('Energy transfer during deceleration');
xlabel('\tau_{CVT} ');
ylabel('\DeltaE_{Kers} [J]');
legend({num2str(tauk')})

figure()
colormap('jet')
[X,Y] = meshgrid(tauCVT, tauk);
Z = T0*Jf/M*(Jeqv/(Jeqv+Jeqk))^2*(1./(tauf*taugb*Rw*X.*Y)).^2;
surf(X,Y,Z)
xlabel('\tau_{CVT}')
ylabel('\tau_k')
title('Energy transfer with variables \tau_{CVT}, \tau_k');
zlabel('\DeltaE_{Kers} [J]')


figure()
for i=1:length(tauk)
    omegaf = sqrt(2*E(i,:)/Jf);
    plot(tauCVT, omegaf*9.55);
    %text(v(100)*3.6, omegaf(100)*9.55, num2str(gamma(i)));
    hold on;
end
plot(tauCVT, 64.5*10^3*ones(1, 100), 'black--');
%ylim([0, 8.5*10^4]);
xlim([taumin, taumax]);
title('Velocity comparison during deceleration');
xlabel('\tau_{CVT}');
ylabel('\omega_f [rpm]');
legend({num2str(tauk')});