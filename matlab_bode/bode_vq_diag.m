%%% Bode plots for (rho,q)
%% %% distributed, lambda2 > 0
clear; clc; close all;

% define parameters
L = 100; % section length, [m]
tau = 15; % relaxation time, [s]
rho0 = 0.01; % rho*, linearization point 
minLogFreq = -4;
maxLogFreq = 2;

% Greenshields Hamiltonian 
rhomax = 0.1; % jam density, [veh/m]
qmax = 1300/3600; % max flow, [veh/s]

q = @(rho) 4*qmax*rho*(rhomax - rho)/(rhomax^2);
qprime = @(rho) 4*qmax/rhomax - 8*qmax*rho/rhomax^2;

% Triangular Hamiltonian
% rhomax = 0.1274; % jam density, [veh/m]
% qmax = 2040/3600; % max flow, [veh/s]
% rhoc = 0.0186; % critical density, [veh/m]
% 
% q = @(rho) (rho<=rhoc)*(qmax*rho/rhoc) + (rho>rhoc)*(qmax*(rho-rhomax)/(rhoc - rhomax));
% qprime = @(rho) (rho<=rhoc)*(qmax/rhoc) + (rho>rhoc)*qmax/(rhoc - rhomax);

q0 = q(rho0);
lambda1 = q0/rho0 ; % lambda1 = v* = q(rho*)/rho*
lambda2 = qprime(rho0); % lambda2 = v* + rho* V'(rho*) = q'(rho*)

alpha = - lambda2 / (tau * (lambda1 - lambda2))
cutoff = 2 * pi * lambda1 * tau * alpha / L

w = logspace(minLogFreq,maxLogFreq,1000); % frequency points
s = 2*pi*1i*w;
X = linspace(0, L, 1000);
X(end) = []; % x=L causes artifact

phi11 = NaN(length(X),length(s));
phi21 = phi11;
phi22 = phi11;

for i = 1:length(X);
    x = X(i);
    phi11(i,:) = exp((-x/lambda1)*(s+1/tau));
    phi21(i,:) = lambda1*(exp((-x/lambda1)*(s+1/tau)) - exp(-x*s/lambda2))./...
        (lambda2 - s*tau*(lambda1-lambda2));
    phi22(i,:) = exp(-x*s/lambda2);
end
dB11 = 20*log10(abs(phi11));
phase11 = (180/pi)*unwrap(angle(phi11));
dB21 = 20*log10(abs(phi21));
phase21 = (180/pi)*unwrap(angle(phi21));
dB22 = 20*log10(abs(phi22));
phase22 = (180/pi)*unwrap(angle(phi22));
%% Spatial Bode, lambda2 > 0
fig1 = figure(1);
set(fig1,'defaulttextinterpreter','latex');
surf(w,X,dB11,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
zlabel('gain [dB]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
title('Bode plot for $\phi_{11}(x,s)$')
view([1 -2 1])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
print(fig1,'-dpdf','distr_phi_11')

fig2 = figure(2);
set(fig2,'defaulttextinterpreter','latex');
surf(w,X,dB21,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
zlabel('gain [dB]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
title('Bode plot for $\phi_{21}(x,s)$')
view([1 -1.8 2])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
print(fig2,'-dpdf','distr_phi_21')

%% distributed, lambda2 < 0
clear; close all;

% define parameters
L = 100; % section length, [m]
tau = 15; % relaxation time, [s]
rho0 = 0.08; % rho*, linearization point
minLogFreq = -4;
maxLogFreq = 2;

% Greenshields Hamiltonian 
rhomax = 0.1; % jam density, [veh/m]
qmax = 1300/3600; % max flow, [veh/s]

q = @(rho) 4*qmax*rho*(rhomax - rho)/(rhomax^2);
qprime = @(rho) 4*qmax/rhomax - 8*qmax*rho/rhomax^2;

% Triangular Hamiltonian
% rhomax = 0.1274; % jam density, [veh/m]
% qmax = 2040/3600; % max flow, [veh/s]
% rhoc = 0.0186; % critical density, [veh/m]
% 
% q = @(rho) (rho<=rhoc)*(qmax*rho/rhoc) + (rho>rhoc)*(qmax*(rho-rhomax)/(rhoc - rhomax));
% qprime = @(rho) (rho<=rhoc)*(qmax/rhoc) + (rho>rhoc)*qmax/(rhoc - rhomax);

q0 = q(rho0);
lambda1 = q0/rho0 ; % lambda1 = v* = q(rho*)/rho*
lambda2 = qprime(rho0); % lambda2 = v* + rho* V'(rho*) = q'(rho*)

alpha = - lambda2 / (tau * (lambda1 - lambda2))
cutoff = 2 * pi * lambda1 * tau * alpha / L

w = logspace(minLogFreq,maxLogFreq,800); % frequency points
s = 1i*w;
X = 0:.1:L;
X(end) = []; % x=L causes artifact

gamma11 = NaN(length(X),length(s));
gamma21 = gamma11;
gamma22 = gamma11;

for i = 1:length(X);
    x = X(i);
    gamma11(i,:) = exp((-x/lambda1)*(s+1/tau));
    gamma21(i,:) = lambda1*(exp((-x/lambda1)*(s+1/tau)) - ...
        exp(-L/(tau*lambda1) - (x - L*(lambda1 - lambda2)/lambda1)*s/lambda2))./...
        (lambda2 - s*tau*(lambda1-lambda2));
    gamma22(i,:) = exp(-(x-L)*s/lambda2);
end
dB11 = 20*log10(abs(gamma11));
phase11 = (180/pi)*unwrap(angle(gamma11));
dB21 = 20*log10(abs(gamma21));
phase21 = (180/pi)*unwrap(angle(gamma21));
dB22 = 20*log10(abs(gamma22));
phase22 = (180/pi)*unwrap(angle(gamma22));
%% Spatial Bode, lambda2 < 0
fig1 = figure(1);
set(fig1,'defaulttextinterpreter','latex');
surf(w,X,dB11,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
title('Bode plot for $\gamma_{11}(x,s)$')
view([1 -2 1])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
print(fig1,'-dpdf','distr_gamma_11')

fig2 = figure(2);
set(fig2,'defaulttextinterpreter','latex');
surf(w,X,dB21,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
title('Bode plot for $\gamma_{21}(x,s)$')
view([1 -1.8 2])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
print(fig2,'-dpdf','distr_gamma_21')