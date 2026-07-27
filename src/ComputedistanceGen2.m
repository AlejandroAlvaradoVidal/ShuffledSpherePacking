function [d] = ComputedistanceGen2(G,NF,R)
%G=G(1:NF,1:NF,:);
[r,c,z] = find(G);
X = [r c z]'; % generator matrix
AllPoints = X';% your points
%R = 1; % depends on your data
for idx=1:length(AllPoints)
Distances = sqrt( sum( (AllPoints-AllPoints(idx,:)).^2 ,2) );
Dis=sort(Distances);
Ninside   = length( find(Dis<=R) );
end
d = Dis(Ninside);
end