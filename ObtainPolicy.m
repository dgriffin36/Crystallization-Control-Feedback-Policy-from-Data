function [Policy,Grid]=ObtainPolicy(Xtr,Utr,dXtr,dt,Dt,xTarget,N,rho,gamma,Grid)             
% ObtainPolicy gives a state-feedback control policy for controlling
% batch cooling crystallization according to the method described in
% "Data-Driven Modeling and Dynamic Programming Applied to Batch Cooling
% Crystallization" by D. J. Griffin, M. A. Grover, Y. Kawajiri, and R. W.
% Rousseau.
%
% To run this function, the m-files in the subfunctions folder must be in 
% the same path. In addition, CVX (Software for Disciplined Convex 
% Programming [1],[2]) is required. This must be installed and in the 
% appropriate path: http://cvxr.com/cvx/download/. Note: for commercial 
% use with non-free solvers, such as MATLAB, please obtain the appropriate 
% license (http://cvxr.com/cvx/licensing/).
%
%-----------------------------------Inputs---------------------------------
% There are a number of inputs. The inputs specify: the training 
% data, the length of the time steps, the target position in mass-count 
% space, the batch run time, adjustable parameters in the optimization 
% formulation, and the space-input discretization to use.
%   
% REQUIRED
%   Xtr     - Contains training data 'positions'.
%   Utr     - Contains training data inputs.
%   dXtr    - Contains training data 'movements'. 
%   dt      - The measurement time interval in minutes.
%   Dt      - The time interval over which the input is held constant.
%   xTarget - The target position [count, mass].
%   N       - The number of time intervals in the control run.
%
% OPTIONAL
%   rho     - Parameter that weights the input-effort cost (default: 25).
%   gamma 	- Parameter that weights the running distance-to-target (5).
%   Grid    - Structure containing the grid spacing for: count (Grid.c), 
%             mass (Grid.m), and supersaturation (Grid.s).
% 
%-------------------------------Note---------------------------------------
% Depending on input data and discretization, the function will require
% substantial computation time (on the order of 30 minutes for the example 
% data). The modeling step takes the longest. Prompts and visuals have been 
% added as progress checks. These should be commented-out to run the 
%f unction unattended.
%
%-----------------------------Bibliography---------------------------------
%[1] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined 
%convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013.
%
%[2] Michael Grant and Stephen Boyd. Graph implementations for nonsmooth 
%convex programs, Recent Advances in Learning and Control (a tribute to M. 
%Vidyasagar), V. Blondel, S. Boyd, and H. Kimura, editors, pages 95-110, 
%Lecture Notes in Control and Information Sciences, Springer, 2008. 
%http://stanford.edu/~boyd/graph_dcp.html.
%
%----------------------------Copyright-------------------------------------
% Copyright 2015, Daniel Griffin.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but without any warranty; without even the implied warranty of
% merchantability or fitness for a particular purpose. See the GNU General 
% Public License for more details <http://www.gnu.org/licenses/>.


%% check on necessary inputs
if nargin<7
    disp('not enough inputs')
    return
end

%% check for optimization parameters and set defaults if unspecified
if nargin < 8
    rho = 25; 
end

if nargin < 9
    gamma = 5;
end

%% check for specification of the discretization scheme
if nargin < 10
    % ---- calculate default discretization ----------------
    
    % state-space bounds
    countUL = max(Xtr(:,1));
    massUL = max(Xtr(:,2));

    % supersaturation-only model of the dynamics to gauge movement
    k=1;
    for i=1:length(Utr)
        if Xtr(i,1) > 25 && Xtr(i,2) > 1
            X(k,:)=[Xtr(i,:) Utr(i)];
            Y(k,:)=dXtr(i,:);
            k=k+1;
        end
    end
    
    minSup=min(X(:,3)); 
    maxSup=max(X(:,3));
    sGrid=[minSup:abs(minSup)/3:0, 0+maxSup/6:maxSup/6:maxSup];
    ns = length(sGrid);

    % linear regression
    Bc=[X(:,3) X(:,3).^2 X(:,3).^3 X(:,3).^4 X(:,3).^5 X(:,3).^6]\Y(:,1);
    Bm=[X(:,3) X(:,3).^2 X(:,3).^3 X(:,3).^4 X(:,3).^5 X(:,3).^6]\Y(:,2);

    % set grid spacing based on movement according to the sup-only model.
    sn1=sGrid(2);
    sp1=sGrid(6);
    dcn1=[sn1 sn1^2 sn1^3 sn1^4 sn1^5 sn1^6]*Bc;
    dcp1=[sp1 sp1^2 sp1^3 sp1^4 sp1^5 sp1^6]*Bc;
    dmn1=[sn1 sn1^2 sn1^3 sn1^4 sn1^5 sn1^6]*Bm;
    dmp1=[sp1 sp1^2 sp1^3 sp1^4 sp1^5 sp1^6]*Bm;
    cSpacing=min([abs(dcn1)*(Dt/dt)/2;dcp1*(Dt/dt)/2]);
    mSpacing=min([abs(dmn1)*(Dt/dt)/2;dmp1*(Dt/dt)/2]);
    cGrid=0:cSpacing:countUL;
    mGrid=0:mSpacing:massUL;
    nm=length(mGrid);
    nc=length(cGrid);
    Grid.c=cGrid;
    Grid.m=mGrid;
    Grid.s=sGrid;
    scale=cSpacing/mSpacing;
    Grid.scale=scale;
else
    % unpack grid scheme if specified
    cGrid=Grid.c;
    mGrid=Grid.m;
    sGrid=Grid.s;
    scale=Grid.scale;
    nc=length(cGrid);
    nm=length(mGrid);
    ns=length(sGrid);
    countUL=max(cGrid);
    massUL=max(mGrid);
end

%% -------------- graph for visual check -----------------------------------
% visualize the spacing
figure(1)
subplot(1,2,1)
for j=1:length(cGrid)
    plot(cGrid(j)*ones(1,length(mGrid)),mGrid,'.b');hold on;
end
xlabel('chord count')
ylabel('crystal mass [g]')
title('state-space grid')
axis([0 countUL 0 massUL])
subplot(1,2,2)
plot(1:length(sGrid),sGrid,'ob')
xlabel('input level index')
ylabel('supersaturation')
title('input discretization')
% ------------end graph for visual check ----------------------------------

%% comment-out the following prompt to run the function unattended --------
prompt = 'Grid scheme set. Proceed to bandwidth calculations? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end

if str == 'N'
    return
end
%-------------------------end prompt--------------------------------------%

%% obtain bandwidth for each grid point (bandwidth used in the learning algorithm weight function)

kappa=zeros(nc,nm);

for i = 1:nm
    for j = 1:nc
        kappa(j,i)=findKappa(cGrid(j),mGrid(i),scale,[Xtr, Utr],0.25); 
        % findKappa determines the bandwidth according to the data density
        % (used in the non-parametric learning algorithm)
    end
end

figure(2)
meshc(cGrid,mGrid,kappa');
xlabel('chord count','FontSize',14,'Interpreter','latex')
ylabel('crystal mass [g]','FontSize',14,'Interpreter','latex')
zlabel('bandwidth $\kappa$','FontSize',14,'Interpreter','latex')

%% comment-out the following prompt to run the function unattended -------
prompt = 'Bandwidth calculated. Proceed to model development? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end

if str == 'N'
    return
end
%-------------------------end prompt--------------------------------------%
%% obtain dynamic model

if length(Xtr(:,1))<1000
%% --comment-out the following prompt to run the function unattended -----
    prompt = 'Small Training Set. Proceed with model development? Y/N [Y]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'Y';
    end

    if str == 'N'
        return
    end  
% -------------------------end prompt-------------------------------------
end
q=0;
Bc=[];Bm=[];
for i=1:nm
    if i/nm>=.25 && q<.5
        disp('model 25% complete')
        q=1;
    elseif i/nm>=.5 && q<1.5
        disp('model 50% complete')
        q=2;
    elseif i/nm>=.75 && q<2.5
        disp('model 75% complete')
        q=3;
    end
        
    for j=1:nc
        s=(j-1)*nm+i;
        if length(Xtr(:,1))>1000
            [Bc(:,s),Bm(:,s)]=LocalModelCoefficients(Xtr,Utr,dXtr,cGrid(j),mGrid(i),kappa(j,i),scale);
            %this function utilizes CVX and is a constrained optimization
        else
            [Bc(:,s),Bm(:,s)]=LocalModelCoefficientsUnconstrained(Xtr,Utr,dXtr,cGrid(j),mGrid(i),kappa(j,i),scale);
            %this function is an unconstrained optimization and may lead to
            %an inaccurate model
        end
    end
end

%% -------------- graph for visual check -----------------------------------
% quiver plots of the model
Dc=zeros(nc,nm,ns);
Dm=zeros(nc,nm,ns);
for i = 1:nm
    for j = 1:nc
        for h=1:ns
            if length(Xtr(:,1))>1000
                sup=[sGrid(h),sGrid(h)^2,sGrid(h)^3,sGrid(h)^4,sGrid(h)^5,sGrid(h)^6];
                Dc(j,i,h)=sup*Bc(:,(j-1)*nm+i);
                Dm(j,i,h)=sup*Bm(:,(j-1)*nm+i);
            else
                sup=[sGrid(h),sGrid(h)^2,sGrid(h)^3];
                Dc(j,i,h)=sup*Bc(:,(j-1)*nm+i);
                Dm(j,i,h)=sup*Bm(:,(j-1)*nm+i);
            end
                
        end
    end
end

[x,y]=meshgrid(cGrid,mGrid);
%
figure1 = figure;

% Create subplot
%_______ 1 ______ 4
subplot1 = subplot(2,3,1,'Parent',figure1,...
    'YTick',[0 4 8 12 16],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 4 8 12 16],...
    'TickLabelInterpreter','latex',...
    'FontSize',12,...
    'FontName','Calibri');

% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot1,[0 16]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(subplot1,[0 18]);
box(subplot1,'off');
hold(subplot1,'on');

% Create quiver
quiver(x/25, y,(Dc(:,:,5)/25)',Dm(:,:,5)','LineWidth',1,'Color',[0 0 0],'AutoScale','off');
title(['u = ' num2str(sGrid(5))],'FontSize',14,'Interpreter','latex');
xlabel('chord count ($x_1$)','Interpreter','latex');
ylabel('crystal mass ($x_2$) [g]','Interpreter','latex');

%_______ 2 ______ 5
subplot1 = subplot(2,3,2,'Parent',figure1,...
    'YTick',[0 4 8 12 16],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 4 8 12 16],...
    'TickLabelInterpreter','latex',...
    'FontSize',12,...
    'FontName','Calibri');
% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot1,[0 16]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(subplot1,[0 18]);
box(subplot1,'off');
hold(subplot1,'on');

% Create quiver
quiver(x/25, y,(Dc(:,:,6)/25)',Dm(:,:,6)','LineWidth',1,'Color',[0 0 0],'AutoScale','off');
title(['u = ' num2str(sGrid(6))],'FontSize',14,'Interpreter','latex');
xlabel('chord count ($x_1$)','Interpreter','latex');
ylabel('crystal mass ($x_2$) [g]','Interpreter','latex');

%_______ 3 ______ 6
subplot1 = subplot(2,3,3,'Parent',figure1,...
    'YTick',[0 4 8 12 16],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 4 8 12 16],...
    'TickLabelInterpreter','latex',...
    'FontSize',12,...
    'FontName','Calibri');
% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot1,[0 16]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(subplot1,[0 18]);
box(subplot1,'off');
hold(subplot1,'on');

% Create quiver
quiver(x/25, y,(Dc(:,:,7)/25)',Dm(:,:,7)','LineWidth',1,'Color',[0 0 0],'AutoScale','off');
title(['u = ' num2str(sGrid(7))],'FontSize',14,'Interpreter','latex');
xlabel('chord count ($x_1$)','Interpreter','latex');
ylabel('crystal mass ($x_2$) [g]','Interpreter','latex');

%_______ 4 ______ 1
subplot1 = subplot(2,3,4,'Parent',figure1,...
    'YTick',[0 4 8 12 16],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 4 8 12 16],...
    'TickLabelInterpreter','latex',...
    'FontSize',12,...
    'FontName','Calibri');
% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot1,[0 16]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(subplot1,[0 18]);
box(subplot1,'off');
hold(subplot1,'on');

% Create quiver
quiver(x/25, y,(Dc(:,:,1)/25)',Dm(:,:,1)','LineWidth',1,'Color',[0 0 0],'AutoScale','off');
title(['u = ' num2str(sGrid(1))],'FontSize',14,'Interpreter','latex');
xlabel('chord count ($x_1$)','Interpreter','latex');
ylabel('crystal mass ($x_2$) [g]','Interpreter','latex');

%_______ 5 ______ 2
subplot1 = subplot(2,3,5,'Parent',figure1,...
    'YTick',[0 4 8 12 16],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 4 8 12 16],...
    'TickLabelInterpreter','latex',...
    'FontSize',12,...
    'FontName','Calibri');
% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot1,[0 16]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(subplot1,[0 18]);
box(subplot1,'off');
hold(subplot1,'on');

% Create quiver
quiver(x/25, y,(Dc(:,:,2)/25)',Dm(:,:,2)','LineWidth',1,'Color',[0 0 0],'AutoScale','off');
title(['u = ' num2str(sGrid(2))],'FontSize',14,'Interpreter','latex');
xlabel('chord count ($x_1$)','Interpreter','latex');
ylabel('crystal mass ($x_2$) [g]','Interpreter','latex');

%_______ 6 ______ 3
subplot1 = subplot(2,3,6,'Parent',figure1,...
    'YTick',[0 4 8 12 16],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 4 8 12 16],...
    'TickLabelInterpreter','latex',...
    'FontSize',12,...
    'FontName','Calibri');
% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot1,[0 16]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(subplot1,[0 18]);
box(subplot1,'off');
hold(subplot1,'on');

% Create quiver
quiver(x/25, y,(Dc(:,:,3)/25)',Dm(:,:,3)','LineWidth',1,'Color',[0 0 0],'AutoScale','off');
title(['u = ' num2str(sGrid(3))],'FontSize',14,'Interpreter','latex');
xlabel('chord count ($x_1$)','Interpreter','latex');
ylabel('crystal mass ($x_2$) [g]','Interpreter','latex');
hold off;

% ------------end graph for visual check ----------------------------------

%% comment-out the following prompt to run the function unattended --------
prompt = 'Model Constructed. Proceed to policy development? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end

if str == 'N'
    return
end
%--------------------------end prompt--------------------------------------
%% convert to cell-to-cell mapping

% obtain the 'next-state' matrix
NS=zeros(nm*nc,ns);
delta_t=Dt/dt;
if length(Xtr(:,1))>1000
    model=1;
else
    model=2;
end

for i=1:nm
    for j=1:nc
        for k=1:ns;
            NS((j-1)*nm+i,k)=IndexNextState([cGrid(j);mGrid(i)],sGrid(k),(j-1)*nm+i,delta_t,Bc,Bm,cGrid,mGrid,model);
            % this function calculates a matrix that specifies the next state
            % according to the dynamic model obtained previously under each
            % different input.
        end
    end
end
%% Dynamic Programming to obtain policy
% once the next-state matrix is obtained, use dynamic programming to find
% the policy
Policy=DynamicProgramming(N,xTarget,gamma,rho,NS,cGrid,mGrid,sGrid,scale); 

%% -------------- graph for visual check -----------------------------------
% view policy as time-variant color map
sTar=xTarget;
Asup=Policy;
t=N-1;
acur=reshape(Asup(:,t),nm,nc);
    fig = figure('Colormap',...
    [0.208100005984306 0.190400004386902 0.578495264053345;...
    0.208100005984306 0.21450001001358 0.627790451049805;...
    0.208100005984306 0.238600000739098 0.677085697650909;...
    0.195904761552811 0.264457136392593 0.72790002822876;...
    0.17072856426239 0.291938096284866 0.779247641563416;...
    0.125271424651146 0.324242860078812 0.830271422863007;...
    0.0591333322227001 0.359833329916 0.868333339691162;...
    0.0116952378302813 0.387509524822235 0.881957113742828;...
    0.0059571429155767 0.408614277839661 0.882842838764191;...
    0.0426523797214031 0.464319050312042 0.857011914253235;...
    0.0793476179242134 0.520023822784424 0.831180930137634;...
    0.0788581967353821 0.521970391273499 0.83063542842865;...
    0.0783687829971313 0.523916959762573 0.830089926719666;...
    0.0778793618083 0.525863528251648 0.829544425010681;...
    0.0773899480700493 0.527810096740723 0.828998923301697;...
    0.076900526881218 0.529756605625153 0.828453421592712;...
    0.0764111131429672 0.531703174114227 0.827907919883728;...
    0.0759216919541359 0.533649742603302 0.827362418174744;...
    0.0754322782158852 0.535596311092377 0.826816916465759;...
    0.0749428570270538 0.537542879581451 0.826271414756775;...
    0.472765564918518 0.729555726051331 0.903331756591797;...
    0.87058824300766 0.921568632125854 0.980392158031464;...
    1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;...
    0.850980401039124 0.325490206480026 0.0980392172932625;...
    0.854248404502869 0.352941185235977 0.0816993489861488;...
    0.857516348361969 0.380392163991928 0.0653594806790352;...
    0.860784292221069 0.407843142747879 0.0490196086466312;...
    0.864052295684814 0.43529412150383 0.0326797403395176;...
    0.86732029914856 0.462745100259781 0.0163398701697588;...
    0.87058824300766 0.490196079015732 0;...
    0.896414995193481 0.550495862960815 0.0666619017720222;...
    0.922241747379303 0.610795676708221 0.133323803544044;...
    0.948068499565125 0.671095430850983 0.199985712766647;...
    0.973895251750946 0.731395244598389 0.266647607088089;...
    0.976689457893372 0.735164046287537 0.261066138744354;...
    0.979483604431152 0.73893278837204 0.255484640598297;...
    0.982277810573578 0.742701590061188 0.249903172254562;...
    0.985071957111359 0.746470391750336 0.244321689009666;...
    0.987866163253784 0.750239133834839 0.238740205764771;...
    0.990660309791565 0.754007935523987 0.233158722519875;...
    0.99345451593399 0.757776737213135 0.22757725417614;...
    0.996248662471771 0.761545479297638 0.221995770931244;...
    0.999042868614197 0.765314280986786 0.216414287686348;...
    0.99414050579071 0.784879148006439 0.201217263936996;...
    0.989238083362579 0.804444074630737 0.186020240187645;...
    0.984335720539093 0.824008941650391 0.170823216438293;...
    0.979433298110962 0.843573808670044 0.155626192688942;...
    0.974530935287476 0.863138675689697 0.14042916893959;...
    0.969628572463989 0.882703542709351 0.125232145190239;...
    0.964726150035858 0.902268469333649 0.110035121440887;...
    0.959823787212372 0.921833336353302 0.0948380976915359;...
    0.96435558795929 0.931113958358765 0.162885293364525;...
    0.968887448310852 0.940394639968872 0.230932503938675;...
    0.973419308662415 0.949675261974335 0.298979699611664;...
    0.977951109409332 0.958955883979797 0.367026895284653;...
    0.98248291015625 0.96823650598526 0.435074090957642;...
    0.987014770507812 0.977517127990723 0.503121316432953;...
    0.991546630859375 0.98679780960083 0.571168482303619;...
    0.996078431606293 0.996078431606293 0.639215707778931],...
    'Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',fig,'BoxStyle','full','Layer','top',...
    'CLim',[-0.1331 0.2]);
box(axes1,'on');
hold(axes1,'on');

    contourf(cGrid,mGrid,acur);hold on;plot(sTar(1),sTar(2),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineWidth',1,...
    'LineStyle','none',...
    'Color',[0 0 0]);
    hold on;plot(sTar(1),sTar(2),'MarkerEdgeColor',[0 0 0],'MarkerSize',10,'Marker','+',...
    'LineWidth',1,...
    'LineStyle','none',...
    'Color',[0 0 0]);hold off
    xlabel('chord count')
    ylabel('crystal mass [g]')
    title(['time = ' num2str(t*5) ' min.'])
    colorbar
    % Create textbox
    annotation(fig,'textbox',...
    [0.852365620139969 0.943925233644859 0.124306230559876 0.0433925233644776],...
    'String',{'u'},...
    'FitBoxToText','off',...
    'EdgeColor','none');
% -------------- graph for visual check end -------------------------------

%% comment-out the following prompt to run the function unattended --------
prompt = 'Output Policy? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end

if str == 'N'
    return
end
%--------------------------end prompt--------------------------------------
end
