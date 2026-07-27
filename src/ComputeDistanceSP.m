%% Author: Dr. Nelson Eduardo Diaz Diaz
%% Simulator: Compressive spectral imaging using Sphere Packing
%% Compute sphere packing to generate coded apertures
%% September 10, 2024
function [list,xyz]= ComputeDistanceSP(T)
NF = size(T,3);
X = [];
xyz=[];
for i=1:NF
    tp = T(1:NF,1:NF,i)*i;
    [x,y,z] = find(tp);
    xyz=vertcat(xyz,[x,y,z]);
    temp =  [x,y,z];
    X = [X; temp];
end
%[r,c,z1] = find(G);
%X = [r c z1]'; % generator matrix
X = X';
D = pdist(X');
D1 = squareform(D);
R = zeros(length(X),1);
for i=1:length(X)
     temp = D1(i,:);
     [x2,y2,z2]=find(temp);
     R(i)= min(z2)/2;
end
list=(R(:));
end
