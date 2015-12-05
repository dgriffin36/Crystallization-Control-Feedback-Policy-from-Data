function [sup_sp]=setpointFromPolicy(count,mass,t,Dt,Pol,Grid)
% setpointFromPolicy
% This function gives the supersaturation setpoint for the current
% state according to the policy.
%
%-----------------------------------Inputs---------------------------------
%
%   count   - The current total chord count.
%   mass    - The curernt crystal mass [g].
%   t       - Time from the start of the run [min.].
%   Dt      - Time interval over which the input is held constant [min].
%   Pol     - The state-feedback policy.
%   Grid    - Structure containing the grid spacing for: count (Grid.c), 
%             mass (Grid.m), and supersaturation (Grid.s).

%% 
% find the index for the current state
ind=closestState([count;mass],Grid.c,Grid.m);

% index of time
t_ind=floor(t/Dt);

if t_ind>length(Pol(1,:))-1
    t_ind=length(Pol(1,:))-1;
end

%give the supersaturation according to the policy
sup_sp=Pol(ind,t_ind);


end
