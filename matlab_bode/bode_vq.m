%%% Bode plots for (rho,q)

%% [v(x,s);q(x,s)] = N(s)[v(0,s);q(0,s)], 2D plots

clear; clc; close all;
% define parameters
L = 100; % section length, [m]
tau = 15; % relaxation time, [s]
rho0 = 0.01; % rho*, linearization point

% Greenshields Hamiltonian 
% rhomax = 0.1; % jam density, [veh/m]
% qmax = 1300/3600; % max flow, [veh/s]
% 
% q = @(rho) 4*qmax*rho*(rhomax - rho)/(rhomax^2);
% qprime = @(rho) 4*qmax/rhomax - 8*qmax*rho/rhomax^2;

% Triangular Hamiltonian
rhomax = 0.1274; % jam density, [veh/m]
qmax = 2040/3600; % max flow, [veh/s]
rhoc = 0.0186; % critical density, [veh/m]

q = @(rho) (rho<=rhoc)*(qmax*rho/rhoc) + (rho>rhoc)*(qmax*(rho-rhomax)/(rhoc - rhomax));
qprime = @(rho) (rho<=rhoc)*(qmax/rhoc) + (rho>rhoc)*qmax/(rhoc - rhomax);

q0 = q(rho0);
lambda1 = q0/rho0 ; % lambda1 = v* = q(rho*)/rho*
lambda2 = qprime(rho0); % lambda2 = v* + rho* V'(rho*) = q'(rho*)

w = logspace(-3,-1,800); % frequency points
s = 1i*w;
x = L;

psi11 = (lambda2*exp((-x/lambda1).*(s+1/tau)) - (lambda1 - lambda2)*exp(-x.*s/lambda2)*tau.*s)./...
    (lambda2 - (lambda1 - lambda2)*tau*s);
psi12 = (lambda1 - lambda2)*(exp((-x/lambda1).*(s+1/tau)) - exp(-s.*x/lambda2))./...
    (rho0*(lambda2 - (lambda1 - lambda2)*tau*s));
psi21 = rho0*lambda2*(exp((-x/lambda1).*(s+1/tau)) - exp(-s.*x/lambda2))*tau.*s./...
    (lambda2 - (lambda1 - lambda2)*tau*s);
psi22 = (lambda2*exp(-s.*x/lambda2)- (lambda1 - lambda2)*exp((-x/lambda1).*(s+1/tau))*tau.*s)./...
    (lambda2 - (lambda1 - lambda2)*tau*s);

dB11 = 20*log10(abs(psi11));
phase11 = (180/pi)*unwrap(angle(psi11));
dB12 = 20*log10(abs(psi12));
phase12 = (180/pi)*unwrap(angle(psi12));
dB21 = 20*log10(abs(psi21));
phase21 = (180/pi)*unwrap(angle(psi21));
dB22 = 20*log10(abs(psi22));
phase22 = (180/pi)*unwrap(angle(psi22));

fig1 = figure(1);
subplot(2,2,1)
semilogx(w,dB11)
xlabel('Frequency [Hz]')
ylabel('Gain [dB]')
grid on
xlim([w(1) w(end)])
title('\psi_{11}(L,s)')

subplot(2,2,3)
semilogx(w,phase11)
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
grid on
xlim([w(1) w(end)])

subplot(2,2,2)
semilogx(w,dB12)
xlabel('Frequency [Hz]')
ylabel('Gain [dB]')
grid on
xlim([w(1) w(end)])
title('\psi_{12}(L,s)')

subplot(2,2,4)
semilogx(w,phase12)
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
grid on
xlim([w(1) w(end)])

% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(fig1,'-dpdf','IOv_-3to-1')

fig2 = figure(2)
subplot(2,2,1)
semilogx(w,dB21)
xlabel('Frequency [Hz]')
ylabel('Gain [dB]')
grid on
xlim([w(1) w(end)])
title('\psi_{21}(L,s)')

subplot(2,2,3)
semilogx(w,phase21)
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
grid on
xlim([w(1) w(end)])

subplot(2,2,2)
semilogx(w,dB22)
xlabel('Frequency [Hz]')
ylabel('Gain [dB]')
grid on
xlim([w(1) w(end)])
title('\psi_{22}(L,s)')

subplot(2,2,4)
semilogx(w,phase22)
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
grid on
xlim([w(1) w(end)])


% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(fig2,'-dpdf','IOq_-3to-1')
%% distributed
clear; clc; close all;

% define parameters
L = 100; % section length, [m]
tau = 15; % relaxation time, [s]
rho0 = 0.01; % rho*, linearization point

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

w = logspace(-3,-1,800); % frequency points
s = 1i*w;
X = 0:.1:L;

psi11 = NaN(length(X),length(s));
psi12 = psi11;
psi21 = psi11;
psi22 = psi11;

for i = 1:length(X);
    x = X(i);
    psi11(i,:) = (lambda2*exp((-x/lambda1).*(s+1/tau)) - (lambda1 - lambda2)*exp(-x.*s/lambda2)*tau.*s)./...
    (lambda2 - (lambda1 - lambda2)*tau.*s);
    psi12(i,:) = (lambda1 - lambda2)*(exp((-x/lambda1).*(s+1/tau)) - exp(-s.*x/lambda2))./...
    (rho0*(lambda2 - (lambda1 - lambda2)*tau*s));
    psi21(i,:) = rho0*lambda2*(exp((-x/lambda1).*(s+1/tau)) - exp(-s.*x/lambda2))*tau.*s./...
    (lambda2 - (lambda1 - lambda2)*tau*s);
    psi22(i,:) = (lambda2*exp(-s.*x/lambda2)- (lambda1 - lambda2)*exp((-x/lambda1).*(s+1/tau))*tau.*s)./...
    (lambda2 - (lambda1 - lambda2)*tau*s);
end

dB11 = 20*log10(abs(psi11));
phase11 = (180/pi)*unwrap(angle(psi11));
dB12 = 20*log10(abs(psi12));
phase12 = (180/pi)*unwrap(angle(psi12));
dB21 = 20*log10(abs(psi21));
phase21 = (180/pi)*unwrap(angle(psi21));
dB22 = 20*log10(abs(psi22));
phase22 = (180/pi)*unwrap(angle(psi22));
%%
fig1 = figure(1);
surf(w,X,dB11,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
title('Bode plot for \psi_{11}(x,s)')
view([1 -2 1])
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(fig1,'-dpdf','distr11_-3to-1')

fig2 = figure(2);
surf(w,X,dB12,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
view([1 -1.8 3])
title('Bode plot for \psi_{12}(x,s)')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(fig2,'-dpdf','distr12_-3to-1')

fig3 = figure(3);
surf(w,X,dB21,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
title('Bode plot for \psi_{21}(x,s)')
view([1 -1.8 2])
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(fig3,'-dpdf','distr21_-3to-1')

fig4 = figure(4);
surf(w,X,dB22,'Edgecolor','none')
xlabel('Frequency [Hz]')
ylabel('x [m]')
set(gca,'xscale','log')
set(gca,'Ydir','reverse')
title('Bode plot for \psi_{22}(x,s)')
view([1 -2 1.5])
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(fig4,'-dpdf','distr22_-3to-1')
