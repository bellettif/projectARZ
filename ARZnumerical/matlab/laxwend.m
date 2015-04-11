%% Solving ARZ with Lax-Wendroff
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
% figure(2)
% r = 0:.005:rhomax;
% V = Veq(r);
% plot(r,V)
% title('V_{eq}(\rho)')

% F = @(rho,y) [y + rho.*Veq(rho); y.^2./rho + y.*Veq(rho)]; % scalar inputs
F = @(rho,y) cat(3,y + rho.*Veq(rho),y.^2./rho + y.*Veq(rho));

% mesh
Dx = dx/8;
Dt = Dx/20;
x = 0:Dx:dx*(sp-1);
t = 0:Dt:dt*(ti-1);

% initialize U
U = NaN(length(x),length(t),2);
% IC
U(:,1,1) = interp1q(X',rhomap(:,1),x'); % rho I.C.
U(:,1,2) = U(:,1,1).*(interp1q(X',vmap(:,1),x') - Veq(U(:,1,1))); % y I.C. 
% BC
U(1,:,1) = interp1q(T',rhomap(1,:)',t'); % rho left B.C. 
U(end,:,1) =  interp1q(T',rhomap(end,:)',t'); % rho right B.C.
U(1,:,2) = U(1,:,1).*(interp1q(T',vmap(1,:)',t')' - Veq(U(1,:,1))); % y left B.C.
U(end,:,2) = U(end,:,1).*(interp1q(T',vmap(end,:)',t')' - Veq(U(end,:,1))); % y right B.C. 

tic
mid = 2:length(x)-1;
for n = 2:length(t)
  Uplus = (1/2)*((U(mid+1,n-1,:)+U(mid,n-1,:))-(Dt/Dx)*...
     (F(U(mid+1,n-1,1),U(mid+1,n-1,2))-F(U(mid,n-1,1),U(mid,n-1,2))));
  Uminus = (1/2)*((U(mid,n-1,:)+U(mid-1,n-1,:))-(Dt/Dx)*...
      (F(U(mid,n-1,1),U(mid,n-1,2))-F(U(mid-1,n-1,1),U(mid-1,n-1,2))));
  U(mid,n,:) = U(mid,n-1,:) - (Dt/Dx)*(F(Uplus(:,:,1),Uplus(:,:,2)) ...
      - F(Uminus(:,:,1),Uminus(:,:,2)));
end
toc

rhoLW = U(:,:,1);
vLW = U(:,:,2)./U(:,:,1)+ Veq(U(:,:,1));

%% visualize solution %freezes
tic
hfig = figure(7);
set(hfig,'Position', [100, 100, 1000, 350]);
subplot(1,2,1)
h = pcolor(T,X,rhomap);
set(h, 'EdgeColor', 'none');
title('\rho [veh/m]')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)

subplot(1,2,2)
h = pcolor(t,x,rhoLW);
set(h, 'EdgeColor', 'none');
title('Lax-Wendroff - \rho [veh/m]')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)

% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(hfig,'-dpdf','LWrho')


hfig = figure(8);
set(hfig,'Position', [100, 100, 1000, 350]);
subplot(1,2,1)
h = pcolor(T,X,vmap);
set(h, 'EdgeColor', 'none');
title('v [m/s]')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)

subplot(1,2,2)
h = pcolor(t,x,vLW);
set(h, 'EdgeColor', 'none');
title('Lax-Wendroff - v [m/s]')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(hfig,'-dpdf','LWv')
toc

%% debug 

% tslice = 1;
for tslice = 50:length(T)
    figure(9)
    plot(x,rhoLW(:,round(tslice*dt/Dt)),X,rhomap(:,tslice))
    % plot(x,rhoLW(:,round(tslice*dt/Dt)))
    title(sprintf('\\rho at t = %1.2f s',dt*tslice))
    xlabel('x')
    ylabel('\rho(x)')
    legend('Lax-Wendroff','Data')
    pause(0.2)
end
