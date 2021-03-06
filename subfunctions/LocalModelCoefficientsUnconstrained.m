function [bc, bm]=LocalModelCoefficientsUnconstrained(Xtr,Utr,dXtr,c,m,kappa,scale)
% LocalModelCoefficientsUnconstrained
% This funciton outputs the local model coefficients for the current count
% and mass position given the training data available and the bandwidth.
% The local model does not have physical contraints.

%% 
X=[Xtr,Utr];
Y=dXtr;

w=zeros(length(X(:,1)),1);
for i=1:length(X(:,1))
    w(i)=exp(-(((X(i,1)-c)/scale)^2+(X(i,2)-m)^2)/(kappa));
end

bc=[w.^(1/2).*X(:,3) w.^(1/2).*X(:,3).^2 w.^(1/2).*X(:,3).^3]\(w.^(1/2).*Y(:,1));
bm=[w.^(1/2).*X(:,3) w.^(1/2).*X(:,3).^2 w.^(1/2).*X(:,3).^3]\(w.^(1/2).*Y(:,2));

end