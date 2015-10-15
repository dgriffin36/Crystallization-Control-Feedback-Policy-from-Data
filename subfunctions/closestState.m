function [ind]=closestState(x,cGrid,mGrid)
% closestState 
% Identifies the closest discrete state in the mass-count grid. This is a
% subfunction called by IndexNextState.

nm=length(mGrid);
[~,j]=min(abs(x(1)-cGrid));
[~,i]=min(abs(x(2)-mGrid));
ind=(j-1)*nm+i;

end
