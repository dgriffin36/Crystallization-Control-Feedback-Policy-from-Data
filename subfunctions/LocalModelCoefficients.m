function [bc, bm]=LocalModelCoefficients(Xtr,Utr,dXtr,c,m,kappa,scale)
% LocalModelCoefficients
% Outputs the local model coefficients for the current count
% and mass position given the training data available and the bandwidth.
% This function is a subfunction called by ObtainPolicy.

%-------------------------------Note---------------------------------------
% This function calls CVX: Software for Disciplined Convex 
%Programming [1],[2]. For this function to run, CVX must be installed and 
%in the appropriate path. For commercial use with non-free solvers, such as 
%MATLAB, please obtain the appropriate license 
%(http://cvxr.com/cvx/licensing/).


%% Solve the weighted and constrained least-squares regression 
% this is denoted by \tilde{\mathcal{A}} in the supporting information of 
% the accompanying manuscript

X=[Xtr,Utr];
Y=dXtr;
w=zeros(length(X(:,1)),1);
for i=1:length(X(:,1))
    w(i)=exp(-(((X(i,1)-c)/scale)^2+(X(i,2)-m)^2)/(kappa));
end

%discretized supersaturation domain
sup=(-.15:.001:.5)';
Sup=[sup sup.^2 sup.^3 sup.^4 sup.^5 sup.^6];
dSup=[ones(length(sup),1) 2*sup 3*sup.^2 4*sup.^3 5*sup.^4 6*sup.^5];

% run CVX to solve the optimization problem
cvx_quiet true;
cvx_begin
variable bc(6,1)
minimize(norm([w.^(1/2).*X(:,3) w.^(1/2).*X(:,3).^2 w.^(1/2).*X(:,3).^3 ...
    w.^(1/2).*X(:,3).^4 w.^(1/2).*X(:,3).^5 w.^(1/2).*X(:,3).^6]*bc...
    -w.^(1/2).*Y(:,1)))
subject to 
sup.*(Sup*bc)>=0;
(dSup*bc)>=0;
cvx_end

cvx_begin
variable bm(6,1)
minimize(norm([w.^(1/2).*X(:,3) w.^(1/2).*X(:,3).^2 w.^(1/2).*X(:,3).^3 ...
    w.^(1/2).*X(:,3).^4 w.^(1/2).*X(:,3).^5 w.^(1/2).*X(:,3).^6]*bm...
    -w.^(1/2).*Y(:,2)))
subject to 
sup.*(Sup*bm)>=0;
(dSup*bm)>=0;
cvx_end

end
