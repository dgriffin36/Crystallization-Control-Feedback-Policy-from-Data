function [ind]=closestState(s,cGrid,mGrid)
% closestState 
% Identifies the closest discrete state in the mass-count grid. This is a
% subfunction called by IndexNextState.

nm=length(mGrid);
[~,j]=min(abs(s(1)-cGrid));
[~,i]=min(abs(s(2)-mGrid));
ind=(j-1)*nm+i;

end
