%% Solving ARZ with WENO
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

% mesh
Dx = dx/3;
Dt = Dx/20;
x = 0:Dx:dx*(sp-1);
t = 0:Dt:dt*(ti-1);

% Fundamental diagram (Greenshields)
Veq = @(rho)(4*qmax/rhomax^2)*(rhomax-rho);
Veqpr = -4*qmax/rhomax^2;
% figure(2)
% r = 0:.005:rhomax;
% V = Veq(r);
% plot(r,V)
% title('V_{eq}(\rho)')

% F = @(rho,y) [y + rho.*Veq(rho); y.^2./rho + y.*Veq(rho)]; % scalar inputs
F = @(rho,y) cat(3,y + rho.*Veq(rho),y.^2./rho + y.*Veq(rho));
% Fpr = @(rho,y) cat(3,y + Veq(rho) + rho*Veqpr, -(y./rho).^2 + y.*Veqpr,...
%     rho.*Veq(rho), 2*(y./rho) + Veq(rho));
Fpr1 = @(rho,y) y + Veq(rho) + rho*Veqpr;
Fpr2 = @(rho,y) rho.*Veq(rho);
Fpr3 = @(rho,y) -(y./rho).^2 + y.*Veqpr;
Fpr4 = @(rho,y) 2*(y./rho) + Veq(rho);

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
Aphalf = NaN(2,2,length(x)-1);
Amhalf = NaN(2,2,length(x)-1);
Rplus = NaN(2,2,length(x)-1);
Lplus = NaN(2,2,length(x)-1);
Rminus = NaN(2,2,length(x)-1);
Lminus  = NaN(2,2,length(x)-1);
a_splus = NaN(2,2,length(x)-1);
a_sminus = NaN(2,2,length(x)-1);
mid = 4:length(x)-2;

for n = 1:2%:length(t)-1
    % TVD Runge-Kutta, third-order
    
    % First compute F_{i+1/2}
    % compute A_{i+1/2}
    Unphalf = 0.5*(U(mid+1,n,:)+U(mid,n,:));
    Unphalf1 = Unphalf(:,:,1);
    Unphalf2 = Unphalf(:,:,2);
    Aphalf(1,1,:) = [Fpr1(Unphalf1,Unphalf2) Fpr2(Unphalf1,Unphalf2); 
        Fpr3(Unphalf1,Unphalf2) Fpr4(Unphalf1,Unphalf2)]; % third input is space
    % find left and right eigenvectors of A_{i+1/2}
    for i = 1:length(x)-1
        [Rplus(:,:,i),a_splus(:,:,i),Lplus(:,:,i)] = eig(Aphalf(:,:,i));
    end
    a_splus1 = max(abs(a_splus(1,1,:)));
    a_splus2 = max(abs(a_splus(2,2,:)));
    % compute F_{i+1/2}
    Fphalf = (1/12)*(-F(U(mid-1,n,1),U(mid-1,n,2)) ...
        + 7*F(U(mid,n,1),U(mid,n,2)) ...
        + 7*F(U(mid+1,n,1),U(mid+1,n,2))...
        - F(U(mid+2,n,1),U(mid+2,n,2))) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,1,1,mid,n,a_splus1,Lplus,U) ...
        + phiN(5/2,3/2,1/2,-1/2,1,-1,mid,n,a_splus1,Lplus,U))*Rplus(1,1,mid) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,2,1,mid,n,a_splus2,Lplus,U) ...
        + phiN(5/2,3/2,1/2,-1/2,2,-1,mid,n,a_splus2,Lplus,U))*Rplus(2,2,mid);
    
    % Again for F_{i-1/2}
    % compute A_{i-1/2}
    Unmhalf = 0.5*(U(mid,n,:)+U(mid-1,n,:));
    Unmhalf1 = Unmhalf(:,:,1);
    Unmhalf2 = Unmhalf(:,:,2);
    Amhalf(1,1,:) = [Fpr1(Unmhalf1,Unmhalf2) Fpr2(Unmhalf1,Unmhalf2); 
        Fpr3(Unmhalf1,Unmhalf2) Fpr4(Unmhalf1,Unmhalf2)];
    % compute left and right eigenvectors of A_{i-1/2}
    for i = 1:length(x)-1
        [Rminus(:,:,i),a_sminus(:,:,i),Lminus(:,:,i)] = eig(Amhalf(:,:,i));
    end
    a_sminus1 = max(abs(a_sminus(1,1,:)));
    a_sminus2 = max(abs(a_sminus(2,2,:)));
    % compute F_{i-1/2}
    Fmhalf = (1/12)*(-F(U(mid-2,n,1),U(mid-2,n,2)) ...
        + 7*F(U(mid-1,n,1),U(mid-1,n,2)) ...
        + 7*F(U(mid,n,1),U(mid,n,2))...
        - F(U(mid+1,n,1),U(mid+1,n,2))) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,1,1,mid-1,n,a_splus1,Lplus,U) ...
        + phiN(5/2,3/2,1/2,-1/2,1,-1,mid-1,n,a_splus1,Lplus,U))*Rplus(1,1,mid-1) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,2,1,mid-1,n,a_splus2,Lplus,U) ...
        + phiN(5/2,3/2,1/2,-1/2,2,-1,mid-1,n,a_splus2,Lplus,U))*Rplus(2,2,mid-1);
    Un0 = U(:,n,:);
    L0 = -(1/Dx)*(Fphalf - Fmhalf);
    Un1 = Un0 +Dt*L0;
    
    % repeat for step 2
    
    % First compute F_{i+1/2}
    % compute A_{i+1/2}
    Unphalf = 0.5*(Un1(mid+1,n,:)+Un1(mid,n,:));
    Unphalf1 = Unphalf(:,:,1);
    Unphalf2 = Unphalf(:,:,2);
    Aphalf(1,1,:) = [Fpr1(Unphalf1,Unphalf2) Fpr2(Unphalf1,Unphalf2); 
        Fpr3(Unphalf1,Unphalf2) Fpr4(Unphalf1,Unphalf2)]; % third input is space
    % find left and right eigenvectors of A_{i+1/2}
    for i = 1:length(x)-1
        [Rplus(:,:,i),a_splus(:,:,i),Lplus(:,:,i)] = eig(Aphalf(:,:,i));
    end
    a_splus1 = max(abs(a_splus(1,1,:)));
    a_splus2 = max(abs(a_splus(2,2,:)));
    % compute F_{i+1/2}
    Fphalf = (1/12)*(-F(Un1(mid-1,n,1),Un1(mid-1,n,2)) ...
        + 7*F(Un1(mid,n,1),Un1(mid,n,2)) ...
        + 7*F(Un1(mid+1,n,1),Un1(mid+1,n,2))...
        - F(Un1(mid+2,n,1),Un1(mid+2,n,2))) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,1,1,mid,n,a_splus1,Lplus,Un1) ...
        + phiN(5/2,3/2,1/2,-1/2,1,-1,mid,n,a_splus1,Lplus,Un1))*Rplus(1,1,mid) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,2,1,mid,n,a_splus2,Lplus,Un1) ...
        + phiN(5/2,3/2,1/2,-1/2,2,-1,mid,n,a_splus2,Lplus,Un1))*Rplus(2,2,mid);
    
    % Again for F_{i-1/2}
    % compute A_{i-1/2}
    Unmhalf = 0.5*(Un1(mid,n,:)+Un1(mid-1,n,:));
    Unmhalf1 = Unmhalf(:,:,1);
    Unmhalf2 = Unmhalf(:,:,2);
    Amhalf(1,1,:) = [Fpr1(Unmhalf1,Unmhalf2) Fpr2(Unmhalf1,Unmhalf2); 
        Fpr3(Unmhalf1,Unmhalf2) Fpr4(Unmhalf1,Unmhalf2)];
    % compute left and right eigenvectors of A_{i-1/2}
    for i = 1:length(x)-1
        [Rminus(:,:,i),a_sminus(:,:,i),Lminus(:,:,i)] = eig(Amhalf(:,:,i));
    end
    a_sminus1 = max(abs(a_sminus(1,1,:)));
    a_sminus2 = max(abs(a_sminus(2,2,:)));
    % compute F_{i-1/2}
    Fmhalf = (1/12)*(-F(Un1(mid-2,n,1),Un1(mid-2,n,2)) ...
        + 7*F(Un1(mid-1,n,1),Un1(mid-1,n,2)) ...
        + 7*F(Un1(mid,n,1),Un1(mid,n,2))...
        - F(Un1(mid+1,n,1),Un1(mid+1,n,2))) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,1,1,mid-1,n,a_splus1,Lplus,Un1) ...
        + phiN(5/2,3/2,1/2,-1/2,1,-1,mid-1,n,a_splus1,Lplus,Un1))*Rplus(1,1,mid-1) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,2,1,mid-1,n,a_splus2,Lplus,Un1) ...
        + phiN(5/2,3/2,1/2,-1/2,2,-1,mid-1,n,a_splus2,Lplus,Un1))*Rplus(2,2,mid-1);
    L1 = -(1/Dx)*(Fphalf - Fmhalf);
    Un2 = Un1 + (Dt/4)*(-3*L0 + L1);
    
    % repeat for step 3
    
    % First compute F_{i+1/2}
    % compute A_{i+1/2}
    Unphalf = 0.5*(Un2(mid+1,n,:)+Un2(mid,n,:));
    Unphalf1 = Unphalf(:,:,1);
    Unphalf2 = Unphalf(:,:,2);
    Aphalf(1,1,:) = [Fpr1(Unphalf1,Unphalf2) Fpr2(Unphalf1,Unphalf2); 
        Fpr3(Unphalf1,Unphalf2) Fpr4(Unphalf1,Unphalf2)]; % third input is space
    % find left and right eigenvectors of A_{i+1/2}
    for i = 1:length(x)-1
        [Rplus(:,:,i),a_splus(:,:,i),Lplus(:,:,i)] = eig(Aphalf(:,:,i));
    end
    a_splus1 = max(abs(a_splus(1,1,:)));
    a_splus2 = max(abs(a_splus(2,2,:)));
    % compute F_{i+1/2}
    Fphalf = (1/12)*(-F(Un2(mid-1,n,1),Un2(mid-1,n,2)) ...
        + 7*F(Un2(mid,n,1),Un2(mid,n,2)) ...
        + 7*F(Un2(mid+1,n,1),Un2(mid+1,n,2))...
        - F(Un2(mid+2,n,1),Un2(mid+2,n,2))) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,1,1,mid,n,a_splus1,Lplus,Un2) ...
        + phiN(5/2,3/2,1/2,-1/2,1,-1,mid,n,a_splus1,Lplus,Un2))*Rplus(1,1,mid) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,2,1,mid,n,a_splus2,Lplus,Un2) ...
        + phiN(5/2,3/2,1/2,-1/2,2,-1,mid,n,a_splus2,Lplus,Un2))*Rplus(2,2,mid);
    
    % Again for F_{i-1/2}
    % compute A_{i-1/2}
    Unmhalf = 0.5*(Un2(mid,n,:)+Un0(mid-1,n,:));
    Unmhalf1 = Unmhalf(:,:,1);
    Unmhalf2 = Unmhalf(:,:,2);
    Amhalf(1,1,:) = [Fpr1(Unmhalf1,Unmhalf2) Fpr2(Unmhalf1,Unmhalf2); 
        Fpr3(Unmhalf1,Unmhalf2) Fpr4(Unmhalf1,Unmhalf2)];
    % compute left and right eigenvectors of A_{i-1/2}
    for i = 1:length(x)-1
        [Rminus(:,:,i),a_sminus(:,:,i),Lminus(:,:,i)] = eig(Amhalf(:,:,i));
    end
    a_sminus1 = max(abs(a_sminus(1,1,:)));
    a_sminus2 = max(abs(a_sminus(2,2,:)));
    % compute F_{i-1/2}
    Fmhalf = (1/12)*(-F(Un2(mid-2,n,1),Un2(mid-2,n,2)) ...
        + 7*F(Un2(mid-1,n,1),Un2(mid-1,n,2)) ...
        + 7*F(Un2(mid,n,1),Un2(mid,n,2))...
        - F(Un2(mid+1,n,1),Un2(mid+1,n,2))) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,1,1,mid-1,n,a_splus1,Lplus,Un2) ...
        + phiN(5/2,3/2,1/2,-1/2,1,-1,mid-1,n,a_splus1,Lplus,Un2))*Rplus(1,1,mid-1) ...
        + (-phiN(-3/2,-1/2,1/2,3/2,2,1,mid-1,n,a_splus2,Lplus,Un2) ...
        + phiN(5/2,3/2,1/2,-1/2,2,-1,mid-1,n,a_splus2,Lplus,Un2))*Rplus(2,2,mid-1);
    L2 = -(1/Dx)*(Fphalf - Fmhalf);
    Un3 = Un3 + (Dt/12)*(-L0 - L1 + 8*L2);
    U(:,n+1,:) = Un3; 
end
