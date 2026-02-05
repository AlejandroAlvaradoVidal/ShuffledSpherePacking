function [d] = ComputedistanceGen(G,N)
%G=G(1:NF,1:NF,:)
[r,c,z] = find(G);
X = [r c z]'; % generator matrix
AllPoints = X';% your points
R = 1; % depends on your data
for idx=1:length(AllPoints)
Distances = sqrt( sum( (AllPoints-AllPoints(idx,:)).^2 ,2) );
Ninside   = length( find(Distances<=R) );
density(idx) = Ninside/(4*pi*R.^3/3);
end
d = min(Ninside);
end