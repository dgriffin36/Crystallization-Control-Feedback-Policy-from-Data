function [Policy]=DynamicProgramming(T,sTar,rho,gamma,NS,cGrid,mGrid,sGrid,scale)
% DynamicProgramming   Applies the dynamic programming algorithm to
% identify a state-feedback policy. Details can be found in Table 1 of
% the accompanying manuscript, "Data-Driven Modeling and Dynamic 
% Programming Applied to Batch Cooling Crystallization." This is  
% subfunction called by ObtainPolicy.

nc=length(cGrid);
nm=length(mGrid);
V=zeros(nm*nc,T);
A=zeros(nm*nc,T-1);
Policy=zeros(nm*nc,T-1);

for i=1:nm
    for j=1:nc
        V((j-1)*nm+i,T)=1*norm([(cGrid(j)-sTar(1))/scale;mGrid(i)-sTar(2)])^2;
    end
end

t=T-1;
while t>=1
    for i=1:nm
        for j=1:nc
            [V((j-1)*nm+i,t),A((j-1)*nm+i,t)]=CalcMinVal([cGrid(j);mGrid(i)],i,j,sTar,V(:,t+1),t,T,rho,gamma,NS,sGrid,mGrid,scale);
            % CalcMinVal indicates the lowest cost-to-go
            Policy((j-1)*nm+i,t)=sGrid(A((j-1)*nm+i,t));
        end
    end
    t=t-1;
end
end  
