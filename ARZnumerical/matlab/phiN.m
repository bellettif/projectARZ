function [ phi ] = phiN(k1,k2,k3,k4,s,pm,mid,n,a_s,L,U)
% phiN is used in WENO approximation at half-grid point
    % time step n, spatial step(s) mid
    % pm = 1 for plus, pm = -1 for minus
    
    %compute \delta F^{s \pm}_{i+k} 
    DF = @(k) 0.5*L(s,s,:)*(F(U(mid+(k-0.5)+1,n,1),U(mid+(k-0.5)+1,n,2)) ...
        - F(U(mid+(k-0.5),n,1),U(mid+(k-0.5),n,2)) ...
        + pm*a_s*(U(mid+(k-0.5)+1,n,:)-U(mid+(k-0.5),n,:)));
    
    % arguments for \phiN
    a = DF(k1);
    b = DF(k2);
    c = DF(k3);
    d = DF(k4);
    
    % smoothness indicators
    IS0 = 13*(a - b).^2 + 3*(a-3*b).^2;
    IS1 = 13*(b - c).^2 + 3*(b + c).^2;
    IS2 = 13*(c - d).^2 + 3*(3*c - d).^2;
    
    ep = 10^-6;
    
    % weights
    a0 = 1/(ep + IS0).^2;
    a1 = 6/(ep + IS1).^2;
    a2 = 3/(ep + IS2).^2;
    
    w0 = a0./(a0 + a1 + a2);
    w2 = a2./(a0 + a1 + a2);
    
    phi = (1/3)*w0.*(a - 2*b + c) + (1/6)*(w2 - 1/2).*(b - 2*c+d);
end

