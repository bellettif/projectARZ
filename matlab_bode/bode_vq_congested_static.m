%%% Bode plots for (rho,q)

%% [v(x,s);q(x,s)] = N(s)[v(0,s);q(0,s)], 2D plots

clear; clc; close all;
% define parameters
L = 100; % section length, [m]
tau = 15; % relaxation time, [s]
rho0 = 0.08; % rho*, linearization point

% Greenshields Hamiltonian 
rhomax = 0.1; % jam density, [veh/m]
qmax = 1300/3600; % max flow, [veh/s]
q = @(rho) 4*qmax*rho*(rhomax - rho)/(rhomax^2);
qprime = @(rho) 4*qmax/rhomax - 8*qmax*rho/rhomax^2;

% Triangular Hamiltonian
%rhomax = 0.1274; % jam density, [veh/m]
%qmax = 2040/3600; % max flow, [veh/s]
%rhoc = 0.0186; % critical density, [veh/m]

%q = @(rho) (rho<=rhoc)*(qmax*rho/rhoc) + (rho>rhoc)*(qmax*(rho-rhomax)/(rhoc - rhomax));
%qprime = @(rho) (rho<=rhoc)*(qmax/rhoc) + (rho>rhoc)*qmax/(rhoc - rhomax);

q0 = q(rho0);
lambda1 = q0/rho0 ; % lambda1 = v* = q(rho*)/rho*
lambda2 = qprime(rho0); % lambda2 = v* + rho* V'(rho*) = q'(rho*)
alpha = - lambda2 / (tau * (lambda1 - lambda2));

w = logspace(-3,1,800); % frequency points
s = 1i*w;
x = L;

den = s + alpha * exp(-L/(tau * lambda1 * alpha) * (s + alpha));

theta11 = ...
exp((L - x) .* s / lambda2) ...
.* ...
(s + alpha * exp(-x/(tau * lambda1 * alpha) .* (s + alpha))) ...
./ ...
den;

theta12 = ...
1.0 / (rho0 * tau) ...
* ...
exp(- x .* s / lambda2) ...
.* ...
exp(-L * (s + alpha) / (tau * lambda1 * alpha)) ...
.* ...
(1.0 - exp(-(x - L)/(tau * lambda1 * alpha) .* (s + alpha))) ...
./ ...
den;

theta21 = ...
rho0 * s * alpha * tau ...
.* ...
exp((L - x) .* s / lambda2) ...
.* ...
(1.0 - exp(-x/(tau * lambda1 * alpha) .* (s + alpha))) ...
./ ...
den;

theta22 = ...
exp(- x .* s / lambda2) ...
.* ...
exp(-L * (s + alpha) / (tau * lambda1 * alpha)) ...
.* ...
(alpha + s .* exp(-(x - L)/(tau * lambda1 * alpha) .* (s + alpha))) ...
./ ...
den;

dB11 = 20*log10(abs(theta11));
phase11 = (180/pi)*unwrap(angle(theta11));
dB12 = 20*log10(abs(theta12));
phase12 = (180/pi)*unwrap(angle(theta12));
dB21 = 20*log10(abs(theta21));
phase21 = (180/pi)*unwrap(angle(theta21));
dB22 = 20*log10(abs(theta22));
phase22 = (180/pi)*unwrap(angle(theta22));

fig1 = figure(1);
set(fig1,'defaulttextinterpreter','latex');
subplot(2,2,1)
semilogx(w,dB11)
xlabel('Frequency [Hz]')
ylabel('Gain [dB]')
grid on
xlim([w(1) w(end)])
title('$\theta_{11}(L,s)$')

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
title('$\theta_{12}(L,s)$')

subplot(2,2,4)
semilogx(w,phase12)
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
grid on
xlim([w(1) w(end)])

set(findall(gcf,'-property','FontSize'),'FontSize',14)
print(fig1,'-dpdf','IO_theta_11_12')

fig2 = figure(2);
set(fig2,'defaulttextinterpreter','latex');
subplot(2,2,1)
semilogx(w,dB21)
xlabel('Frequency [Hz]')
ylabel('Gain [dB]')
grid on
xlim([w(1) w(end)])
title('$\theta_{21}(L,s)$')

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
title('$\theta_{12}(L,s)$')

subplot(2,2,4)
semilogx(w,phase22)
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
grid on
xlim([w(1) w(end)])

set(findall(gcf,'-property','FontSize'),'FontSize',14)
print(fig1,'-dpdf','IO_theta_21_22')