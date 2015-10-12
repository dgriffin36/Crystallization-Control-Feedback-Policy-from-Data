function kappa=findKappa(c,m,Xtrain,fr)
% findKappa Identifies the appropriate bandwidth for local regression.
%           This is a subfunction called by ObtainPolicy. The output is
%           prerequisite to obtaining the dynamic model.           

N=length(Xtrain(:,1));
scale=25;
kappa=500;
U=50000;
L=0;
e=1;

while e>.05
    frac=sum(exp(-(((Xtrain(:,1)-c)/scale).^2+(Xtrain(:,2)-m).^2)/(kappa)))/N;
    if frac>fr+.05
        U=kappa;
        kappa=(U+L)/2;
    elseif frac < fr-.05
        L=kappa;
        kappa=(U+L)/2;
    end
    e=abs(frac-fr);  
end
