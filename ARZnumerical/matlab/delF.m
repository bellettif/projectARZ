function [ DF ] = delF(k,s,pm,n,mid,a_s,L,U)
%delF evaluates \delta F^{s \pm}_{i+k} 
%   at time step n for spatial step(s) mid
%   pm = 1 for plus, pm = -1 for minus
    m = k-0.5;
    DF = 0.5*L(s,s,:)*(F(U(mid+m+1,n,1),U(mid+m+1,n,2)) ...
        - F(U(mid+m,n,1),U(mid+m,n,2)) ...
        + pm*a_s*(U(mid+m+1,n,:)-U(mid+m,n,:)));

end

