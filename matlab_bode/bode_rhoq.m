%%% Bode plots for (rho,q)

%% [rho(0,s);rho(L,s)] = P(s)[q(0,s);q(L,s)]

clear; clc; 
% define parameters
L = ; % section length
alpha = ; % alpha = lambda1 = v*
beta = ; % beta = -lambda2 = v* + rho* V'(rho*)
tau = ; % relaxation time
sigma = -1/tau; 
delta = -beta/tau;

w = logspace(-5,-1,500); % frequency points
s = 1i*w;

d = (s*(alpha - beta) + delta).^2 - 4*alpha*beta*s.*(s-sigma);
eta1 = (s*(alpha - beta) + delta - sqrt(d))/(2*alpha*beta);
eta2 = (s*(alpha - beta) + delta + sqrt(d))/(2*alpha*beta);

p11 = (eta2.*exp(eta1*L) - eta1.*exp(eta2*L))./(s.*(exp(eta2*L) - exp(eta1*L)));
p12 = (eta1 - eta2)./(s.*(exp(eta2*L) - exp(eta1*L)));
p21 = (eta2 - eta1).*exp((eta1+eta2)*L)./(s.*(exp(eta2*L) - exp(eta1*L)));
p22 = (eta1.*exp(eta1*L) - eta2.*exp(eta2*L))./(s.*(exp(eta2*L) - exp(eta1*L)));


dB11 = 20*log10(abs(p11));
phase11 = (180/pi)*unwrap(angle(p11));
dB12 = 20*log10(abs(p12));
phase12 = (180/pi)*unwrap(angle(p12));
dB21 = 20*log10(abs(p21));
phase21 = (180/pi)*unwrap(angle(p21));
dB22 = 20*log10(abs(p22));
phase22 = (180/pi)*unwrap(angle(p22));

figure(1)
subplot(2,2,1)
semilogx(w,dB11)
xlabel('Frequency [rad/s]')
ylabel('Gain [dB]')
grid on
% axis([1e-5 1e-1 -60 20])

subplot(2,2,3)
semilogx(w,phase11)
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
grid on
% axis([1e-5 1e-1 -500 0])

subplot(2,2,2)
semilogx(w,dB12)
xlabel('Frequency [rad/s]')
ylabel('Gain [dB]')
grid on
% axis([1e-5 1e-1 -60 20])

subplot(2,2,4)
semilogx(w,phase12)
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
grid on
% axis([1e-5 1e-1 0 300])

figure(2)
subplot(2,2,1)
semilogx(w,dB21)
xlabel('Frequency [rad/s]')
ylabel('Gain [dB]')
grid on
% axis([1e-5 1e-1 -60 20])

subplot(2,2,3)
semilogx(w,phase21)
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
grid on
% axis([1e-5 1e-1 -500 0])

subplot(2,2,2)
semilogx(w,dB22)
xlabel('Frequency [rad/s]')
ylabel('Gain [dB]')
grid on
% axis([1e-5 1e-1 -60 20])

subplot(2,2,4)
semilogx(w,phase22)
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
grid on
% axis([1e-5 1e-1 0 300])


