function [ind]=IndexNextState(s,a,cur_ind,delta_t,Bc,Bm,cGrid,mGrid)
% IndexNextState    Produces a 'next-state' matrix that specifies the
% cell-to-cell mapping representation of the dynamic model. This is a
% subfunction called by ObtainPolicy that must be run after obtaining the
% dynamic model.

sprev(:,1)=s;
in(1)=cur_ind;
for i=1:delta_t        
    dc=[a a^2 a^3 a^4 a^5 a^6]*Bc(:,in(i));
    dm=[a a^2 a^3 a^4 a^5 a^6]*Bm(:,in(i));
    snew(:,i)=[sprev(1,i)+dc;sprev(2,i)+dm];
    sprev(:,i+1)=snew(:,i);
    in(i+1)=closestState(snew(:,i),cGrid,mGrid);
end
[ind]=closestState(snew(:,delta_t),cGrid,mGrid);
end
