%% Solving ARZ with Lax-Friedrichs

clear; clc; close all;
load('data/parametervariables.mat')
load('data/q_map.mat')
load('data/v_map.mat')
load('data/rho_map.mat')

[sp,ti] = size(qmap);
X = 0:dx:dx*(sp-1);
T = 0:dt:dt*(ti-1);

% hfig = figure(1);
% set(hfig,'Position', [100, 100, 900, 800])
% subplot(2,2,1)
% h = pcolor(T,X,qmap);
% set(h, 'EdgeColor', 'none');
% title('q [veh/s]')
% xlabel('t [s]')
% ylabel('x [m]')
% % colorbar width
% c=colorbar;
% x1=get(gca,'position');
% cpos=get(c,'Position');
% cpos(3)=0.5*cpos(3);
% set(c,'Position',cpos)
% set(gca,'position',x1)
% 
% subplot(2,2,2)
% h = pcolor(T,X,rhomap);
% set(h, 'EdgeColor', 'none');
% title('\rho [veh/m]')
% xlabel('t [s]')
% ylabel('x [m]')
% c=colorbar;
% x1=get(gca,'position');
% cpos=get(c,'Position');
% cpos(3)=0.5*cpos(3);
% set(c,'Position',cpos)
% set(gca,'position',x1)
% 
% subplot(2,2,3)
% h = pcolor(T,X,vmap);
% set(h, 'EdgeColor', 'none');
% title('v [m/s]')
% xlabel('t [s]')
% ylabel('x [m]')
% c=colorbar;
% x1=get(gca,'position');
% cpos=get(c,'Position');
% cpos(3)=0.5*cpos(3);
% set(c,'Position',cpos)
% set(gca,'position',x1)


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
F = @(rho,y, noise) cat(3,y + rho.*Veq(rho),y.^2./rho + y.*Veq(rho) + noise); % gives (time, space, component)


% mesh
Dx = dx;
Dt = Dx/17;
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

% check initalized U
% hfig = figure(3);
% set(hfig,'Position', [100, 100, 1000, 350]);
% subplot(1,2,1)
% plot(X',rhomap(:,1),x',U(:,1,1))
% title('Check \rho(x,0)')
% legend('Data','Interp')
% subplot(1,2,2)
% yinit = rhomap(:,1).*(vmap(:,1)-Veq(rhomap(:,1)));
% plot(X',yinit,x',U(:,1,2))
% title('Check y(x,0)')
% legend('Data','Interp')
% 
% hfig = figure(4);
% set(hfig,'Position', [100, 100, 1000, 350]);
% subplot(1,2,1)
% plot(T',rhomap(1,:),t',U(1,:,1))
% title('Check \rho(0,t)')
% legend('Data','Interp')
% subplot(1,2,2)
% plot(T',rhomap(end,:),t',U(end,:,1))
% title('Check \rho(end,t)')
% legend('Data','Interp')
% 
% hfig = figure(5);
% set(hfig,'Position', [100, 100, 1000, 350]);
% subplot(1,2,1)
% yleft = rhomap(1,:).*(vmap(1,:)-Veq(rhomap(1,:)));
% plot(T',yleft,t',U(1,:,2))
% title('Check y(0,t)')
% legend('Data','Interp')
% subplot(1,2,2)
% yright = rhomap(end,:).*(vmap(end,:)-Veq(rhomap(end,:)));
% plot(T',yright,t',U(end,:,2))
% title('Check y(end,t)')
% legend('Data','Interp')


% noise = @(xidx,tidx) cat(3,zeros(length(x)-2,1),...
%         gaussmf(x(xidx),[std(rhomap(xidx,tidx)) mean(rhomap(xidx,tidx))])');

noise = @(xidx,tidx) 10^56*gaussmf(rhomap(xidx,ceil(tidx*Dt/dt)),...
    [std(qmap(xidx,ceil(tidx*Dt/dt))) mean(qmap(xidx,ceil(tidx*Dt/dt)))]);

tic
mid = 2:length(x)-1;
for n = 1:length(t)-1
    U(mid,n+1,:) = (1/2)*((U(mid+1,n,:)+U(mid-1,n,:))-(Dt/Dx)*...
        (F(U(mid+1,n,1),U(mid+1,n,2),0)-F(U(mid-1,n,1),U(mid-1,n,2),0))); 

    % add noise
%     U(mid,n+1,:) = (1/2)*((U(mid+1,n,:)+U(mid-1,n,:))...
%         -(Dt/Dx)*(F(U(mid+1,n,1),U(mid+1,n,2),0)+ noise(mid+1,ceil(n*Dt/dt))...
%         - F(U(mid-1,n,1),U(mid-1,n,2),0)) - noise(mid-1,ceil(n*Dt/dt)));

%     U(mid,n+1,:) = (1/2)*((U(mid+1,n,:)+cat(3,zeros(length(x)-2,1), noise(mid+1,n)')...
%         + U(mid-1,n,:) + cat(3,zeros(length(x)-2,1), noise(mid-1,n)'))...
%         -(Dt/Dx)*(F(U(mid+1,n,1),U(mid+1,n,2)+noise(mid+1,n)',0)...
%         -F(U(mid-1,n,1),U(mid-1,n,2)+noise(mid+1,n)',0))); 
    

end
toc

rhoLF = U(:,:,1);
vLF = U(:,:,2)./U(:,:,1)+ Veq(U(:,:,1));


% visualize solution
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
h = pcolor(t,x,rhoLF);
set(h, 'EdgeColor', 'none');
title('Lax-Friedrichs - \rho [veh/m]')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)

% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(hfig,'-dpdf','LFrho')

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
h = pcolor(t,x,vLF);
set(h, 'EdgeColor', 'none');
title('Lax-Friedrichs - v [m/s]')
xlabel('t [s]')
ylabel('x [m]')
c=colorbar;
x1=get(gca,'position');
cpos=get(c,'Position');
cpos(3)=0.5*cpos(3);
set(c,'Position',cpos)
set(gca,'position',x1)
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% print(hfig,'-dpdf','LFv')
toc
%% debug 


for tslice = 50:length(T)
    figure(9)
    plot(x,rhoLF(:,round(tslice*dt/Dt)),X,rhomap(:,tslice))
    % plot(x,rhoLF(:,round(tslice*dt/Dt)))
    title(sprintf('\\rho at t = %1.2f s',dt*tslice))
    xlabel('x')
    ylabel('\rho(x)')
    legend('Lax-Friedrichs','Data')
    pause(0.2)
end

%%
figure(10)
for tidx = 1:4:220
    plot(x,rhoLF(:,tidx),X,rhomap(:,2))
    xlabel('x')
    ylabel('\rho(x)')
    legend('Lax-Friedrichs','Data')
    pause(0.05)
end

%% 
clf
figure(11)
% plot(rhomap(mid-1,ceil(n*Dt/dt)),qmap(mid-1,ceil(n*Dt/dt)),'.b')
plot(rhomap(mid-1,ceil(n*Dt/dt)),qmap(mid-1,ceil(n*Dt/dt)),'.b',...
    rhomap(mid-1,ceil(n*Dt/dt)),qmap(mid-1,ceil(n*Dt/dt)) + 10^56*noise(mid-1,n),'.r')
