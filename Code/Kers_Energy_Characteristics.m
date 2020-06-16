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
tauk = 1/3;
%% Clutch parameters
Jc = 0.01; %[m^2*Kg]
%% CVT parameters
taumin = 0.4; 
taumax = 2.3;
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
%% Energy transfer deceleration 

%Generate array with test speeds
v = linspace(1,100);
%v = (250/3.6)*ones(1, 100);

%Generate array with test ratios
%gamma = linspace(taumin/tauk, taumax/tauk, 10);
tauCVT = linspace(taumin, taumax, 10);

%Generate energy curves
E = ones(length(tauCVT),length(v));
E0 = ones(1, length(v));
for i=1:length(tauCVT)
    for j=1:length(v)
        if v(j)<V1max
            taugb = tau1;
        end    
        if v(j)>V1max && v(j)<V2max
            taugb = tau2;
        end    
        if v(j)>V2max && v(j)<V3max
            taugb = tau3;
        end    
        if v(j)>V3max && v(j)<V4max
            taugb = tau4;
        end    
        if v(j)>V4max && v(j)<V5max
            taugb = tau5;
        end    
        if v(j)>V5max && v(j)<V6max
            taugb = tau6;
        end    
        if v(j)>V6max
            taugb = tau7;   
        end
        %taugb = 1;
        %tauf = 1;
        Jeqv = M*(tauf*taugb*Rw/tauCVT(i))^2+Je/(tauCVT(i))^2+4*Jw*(tauf*taugb/tauCVT(i))^2+Jc;
        %Jeqv = M*(Rw*tauf/tauCVT(i))^2+4*Jw*(tauf/tauCVT(i))^2+Jc;
        Jeqk = (Jf/(tauk)^2)+Jc;
        %T0 = 0.5*Jeqv*(tauCVT(i)/tauf*taugb*Rw)^2*v(j)^2;
        T0 = 0.5*Jeqv*(tauCVT(i)/(tauf*Rw))^2*v(j)^2;
        DeltaT = T0*(Jf*Jeqv/(Jeqv+Jeqk)^2)*(1/tauk^2);
        E(i,j) = DeltaT;
        E0(j) = T0;
    end
end

Emax = 593*10^3;

for i=1:length(tauCVT)
    plot(E0, E(i,:)./Emax);
    text(E0(100), E(i, 100)./Emax, num2str(tauCVT(i)));
    hold on;
end
plot(E0, 1*ones(1, 100), 'black--');
title('Energy transfer during deceleration');
xlabel('E_0 [J]');
ylabel('\DeltaE_{Kers} [p.u.]');
%legend({num2str(tauCVT')});


figure()
for i=1:length(tauCVT)
    omegaf = sqrt(2*E(i,:)/Jf);
    plot(v*3.6, omegaf*9.55);
    text(v(100)*3.6, omegaf(100)*9.55, num2str(tauCVT(i)));
    hold on;
end
plot(v*3.6, 64.5*10^3*ones(1, 100), 'black--');
%ylim([0, 8.5*10^4]);
title('Velocity comparison during deceleration');
xlabel('V_0 [Km/h]');
ylabel('\omega_f [rpm]');
%legend({num2str(tauCVT')});

%% Energy transfer Acceleration

% %Generate array with flywheel test speeds
% omegaf = linspace(64500*0.3/9.55, 64500/9.55, 100);
% 
% % Usage limit for the vehicle tang. speed
% % worst case
% limit = (tauk*omegaf)/(taumax)*Rw;
% 
% % Initial velocity definition
% %v = 1*limit;
% v = (250/3.6)*ones(1, 100);
% %Generate array with vehicle test speeds
% % Initial velocity = 0 m/s
% %v = zeros(1,100); %[m/s]
% 
% %Generate array with test ratios
% gamma = linspace(taumin/tauk, taumax/tauk, 10);
% 
% %Generate energy curves
% E = ones(length(gamma),length(v));
% E0 = ones(1, length(v));
% for i=1:length(gamma)
%     for j=1:length(v)
%         if v(j)<V1max
%             taugb = tau1;
%         end    
%         if v(j)>V1max && v(j)<V2max
%             taugb = tau2;
%         end    
%         if v(j)>V2max && v(j)<V3max
%             taugb = tau3;
%         end    
%         if v(j)>V3max && v(j)<V4max
%             taugb = tau4;
%         end    
%         if v(j)>V4max && v(j)<V5max
%             taugb = tau5;
%         end    
%         if v(j)>V5max && v(j)<V6max
%             taugb = tau6;
%         end    
%         if v(j)>V6max
%             taugb = tau7;   
%         end
%         Jeq = 0.5*(Je/(tauf*taugb)+M*Rw^2+4*Jw);
%         phi = 0.5*Jeq*(Jeq^2/(Jc+Jeq)^2-1);
%         psi = Jeq/(Jf)*Jc^2/(Jc+Jeq)^2*(1/gamma(i)^2);
%         rho = Jeq^2*Jc/((Jc+Jeq)^2*gamma(i));
%         E0kers = 0.5*Jf*omegaf(j)^2;
%         DeltaE = phi*(v(j)/Rw)^2+psi*E0kers+rho*(v(j)/Rw)*omegaf(j);
%         E(i,j) = DeltaE;
%         E0(j) = E0kers;
%     end
% end
% 
% for i=1:length(gamma)
%     plot(E0, E(i,:));
%     %text(E0(100), E(i, 100), num2str(gamma(i)));
%     hold on;
% end
% omega_f_min_rpm = 64500*0.3/9.55; %[rad(sec]
% Energy_min = 1/2*Jf*(omega_f_min_rpm)^2;
% xline(Energy_min, 'black--');
% title('Energy transfer during acceleration');
% xlabel('E_{0,Kers} [J]');
% ylabel('\DeltaE_{Vehicle} [J]');
% %ylim([-300 4000])
% legend({num2str(gamma')})
% 
% figure()
% for i=1:length(gamma)
%     DeltaV = sqrt(2*E(i,:).*sign(E(i,:))*Rw^2/Jeq).*sign(E(i,:));
%     plot(omegaf*9.55, DeltaV*3.6);
% %   text(omegaf(100)*9.55, vf(100)*3.6, num2str(gamma(i)));
%     hold on;
% end
% 
% omega_f_min_rpm = 64500*0.3; %[rpm]
% xline(omega_f_min_rpm, 'black--');
% title('Velocity comparison during acceleration');
% xlabel('\omega_{f0} [rpm]');
% ylabel('\DeltaV_{f,Vehicle} [Km/h]');
% %ylim([-5 30]);
% legend({num2str(gamma')});
% 
% % How to plot 3d surfaces
% % [X,Y] = meshgrid(-2:.02:2);                                
% % Z = X.^2+Y.^2;
% % 
% % surf(X,Y,Z)