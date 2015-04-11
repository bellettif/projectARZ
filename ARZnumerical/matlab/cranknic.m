%% Solving LARZ with Crank-Nicolson

clear; clc; close all;
load('data/parametervariables.mat')
load('data/q_map.mat')
load('data/v_map.mat')
load('data/rho_map.mat')

[sp,ti] = size(qmap);
X = 0:dx:dx*(sp-1);
T = 0:dt:dt*(ti-1);

rhomax = max(max(rhomap));
qmax = max(max(qmap));

% Fundamental diagram (Greenshields)
Veq = @(rho)(4*qmax/rhomax^2)*(rhomax-rho);
Veqpr = -4*qmax/rhomax^2;

% mesh
% Dx = dx/3;
% Dt = Dx/45;
Dx = dx/2;
Dt = Dx/40;
x = 0:Dx:dx*(sp-1);
t = 0:Dt:dt*(ti-1);
tic 
% initialize U
U = NaN(2*length(x),length(t));
% (rho, v)
% IC
% U(1:2:end,1) = interp1q(X',rhomap(:,1),x')-rho_star; % rho
% U(1:2:end,1) = smooth(smooth(U(1:2:end,1)));
% U(2:2:end,1) = interp1q(X',vmap(:,1),x')-v_star; % v
% U(2:2:end,1) = smooth(smooth(U(2:2:end,1)));
U(1:2:end,1) = smooth(smooth(interp1q(X',rhomap(:,1),x')-rho_star)); % rho
U(2:2:end,1) = smooth(smooth(interp1q(X',vmap(:,1),x')-v_star)); % v

% BC
% rholeft = interp1q(T',rhomap(1,:)',t')-rho_star; % rho left
% rhoright =  interp1q(T',rhomap(end,:)',t')-rho_star; % rho right
% vleft = interp1q(T',vmap(1,:)',t')-v_star; % y left
% vright = interp1q(T',vmap(end,:)',t')-v_star; % y right

rholeft = (4*cos(t/80)+0.5*cos(t/100)+3*cos(t/40))*rho_star/10; % rho left
rhoright = (4*cos(t/80)+0.2*cos(t/120)+2*cos(t/40))*rho_star/8; % rho right
vleft = (5*cos(t/80)+10*cos(t/40))*v_star/30; % y left
vright = (3*cos(t/80)+8*cos(t/40))*v_star/15; % y right

A = [v_star rho_star;
    0 v_star + rho_star*Veqpr];
R = (Dt/Dx)*A;

A1sup = [zeros(2*length(x)-2,2) kron(eye(length(x)-1),R/4)
    zeros(2,2*length(x))];
A1sub = [zeros(2,2*length(x))
    kron(eye(length(x)-1),-R/4) zeros(2*length(x)-2,2)];
A1 = eye(2*length(x)) + A1sup + A1sub;

A2sup = [zeros(2*length(x)-2,2) kron(eye(length(x)-1),-R/4)
    zeros(2,2*length(x))];
A2sub = [zeros(2,2*length(x))
    kron(eye(length(x)-1),R/4) zeros(2*length(x)-2,2)];
A2 = eye(2*length(x)) + A2sup + A2sub;
% A1inv = inv(A1);

for n = 1:length(t) - 1
    U(:,n+1) = A1\A2*U(:,n);
    % BC's
    U(1,n+1) = rholeft(n+1);
    U(2,n+1) = vleft(n+1);
    U(end-1,n+1) = rhoright(n+1);
    U(end,n+1) = vright(n+1);
end
toc

rhoCN = U(1:2:end,:)+rho_star;
vCN = U(2:2:end,:)+v_star;

%% visualize solution 
tic
hfig = figure(1);
set(hfig,'Position', [100, 100, 1000, 350]);
subplot(1,2,1)
h = pcolor(T,X,rhomap);
set(h, 'EdgeColor', 'none');
title('Data')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
v = caxis;

subplot(1,2,2)
h = pcolor(t,x,rhoCN);
set(h, 'EdgeColor', 'none');
title(sprintf('\\Delta x = %1.4f, \\Delta t = %1.4f',Dx,Dt))
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
caxis(v)
suptitle('Crank-Nicolson - \rho [veh/m]')

% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(hfig,'-dpdf',sprintf('CNrho_dx_%d_dt_%d',3,45)) % numbers refer to ratio: Dx = dx/%d


hfig = figure(2);
set(hfig,'Position', [100, 100, 1000, 350]);
subplot(1,2,1)
h = pcolor(T,X,vmap);
set(h, 'EdgeColor', 'none');
title('Data')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
v = caxis;

subplot(1,2,2)
h = pcolor(t,x,vCN);
set(h, 'EdgeColor', 'none');
title(sprintf('\\Delta x = %1.4f, \\Delta t = %1.4f',Dx,Dt))
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
caxis(v)
suptitle('Crank-Nicolson - v [veh/s]')

% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(hfig,'-dpdf',sprintf('CNrhov_dx_%d_dt_%d',3,45)) % numbers refer to ratio: Dx = dx/%d

toc
%% Debug

figure(10)
for tidx = 1:4:220
    plot(x,rhoCN(:,tidx),X,rhomap(:,1))
    xlabel('x')
    ylabel('\rho(x)')
    legend('Crank-Nicolson','Data')
    title(sprintf('t = %1.2f s',tidx*Dt))
    pause(0.05)
end
%% smooth test
fig = figure(4);
subplot(1,2,1)
plot(x,rhoCN(:,1),'k',X,rhomap(:,1),'r')
legend('Smoothed IC','Data')
title('\rho at t = 0 s')
subplot(1,2,2)
plot(x,vCN(:,1),'k',X,vmap(:,1),'r')
legend('Smoothed IC','Data')
title('v at t = 0 s')

% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% set(fig,'PaperOrientation','landscape');
% set(fig,'PaperUnits','normalized');
% set(fig,'PaperPosition', [0 0.25 1 0.5]);
% print(fig,'-dpdf','smoothICtest')

%% Smooth IC and BC
fig = figure(5);
set(fig,'Position', [100, 100, 1000, 350]);
subplot(1,2,1)
% plot(t,rhoCN(1,:),'k',T,rhomap(1,:),'r',t,rhoCN(end,:),'k',T,rhomap(end,:),'r')
plot(t,rhoCN(1,:),t,rhoCN(end,:))
legend('Left BC','Right BC')
title('\rho')
subplot(1,2,2)
plot(t,vCN(1,:),t,vCN(end,:))
legend('Left BC','Right BC')
title('v')
suptitle('Smooth BCs')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0.25 1 0.5]);
print(fig,'-dpdf','smoothalltest')

hfig = figure(1);
set(hfig,'Position', [100, 100, 1000, 350]);
subplot(1,2,1)
h = pcolor(t,x,rhoCN);
set(h, 'EdgeColor', 'none');
title('Data')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
title('\rho')

subplot(1,2,2)
h = pcolor(t,x,vCN);
set(h, 'EdgeColor', 'none');
title(sprintf('\\Delta x = %1.4f, \\Delta t = %1.4f',Dx,Dt))
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
title('v')
suptitle('Crank-Nicolson (smooth IC & BC) - v [veh/s]')

% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% set(hfig,'PaperOrientation','landscape');
% set(hfig,'PaperUnits','normalized');
% set(hfig,'PaperPosition', [0 0.25 1 0.5]);
% print(hfig,'-dpdf',sprintf('CNsmoothall_dx_%d_dt_%d',2,40)) % numbers refer to ratio: Dx = dx/%d


%%
for tslice = 1:10%length(T)
    figure(9)
    plot(x,rhoCN(:,round(tslice*dt/Dt)),X,rhomap(:,tslice))
    % plot(x,rhoLF(:,round(tslice*dt/Dt)))
    title(sprintf('\\rho at t = %1.2f s',dt*tslice))
    xlabel('x')
    ylabel('\rho(x)')
    legend('Crank-Nicolson','Data')
    pause(2)
end