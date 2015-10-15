function [v,a]=CalcMinVal(x,i,j,xTarget,Vplus,t,N,rho,gamma,NS,sGrid,mGrid,scale)
% CalcMinVal
% Selects the lowest cost-to-go from the current state for each possible
% input. This function is a subfunction called by DynamicProgramming.

ns=length(sGrid);
nm=length(mGrid);
J=zeros(ns,1);

for k=1:ns
    J(k)=(t/N)^gamma*norm([(x(1)-xTarget(1))/scale;x(2)-xTarget(2)])^2+rho*sGrid(k)^2+Vplus(NS((j-1)*nm+i,k));
end
[v,a]=min(J);
end

