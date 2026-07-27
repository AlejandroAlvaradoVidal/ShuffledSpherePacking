function [density,diameter]= ComputeDensityIrregularSP(T)
NF = size(T,3);
is_gpu = isa(T, 'gpuArray');
if is_gpu
    X = gpuArray([]);
else
    X = [];
end

for i=1:NF
    tp = T(1:NF,1:NF,i)*i;
    [x,y,z] = find(tp);
    temp =  [x,y,z];
    X = [X; temp];
end
%[r,c,z1] = find(G);
%X = [r c z1]'; % generator matrix
X = X';
[U,S,V] = svd(X,'econ');
%if(NF <= 48)
Q = round(U*S*V(1:round(length(X)),:)');
%else
%    Q = round(U*S*V(1:round(0.050*length(X)),:)');
%end
D = pdist(Q',"euclidean");
D1 = squareform(D);
if is_gpu
    R = gpuArray.zeros(length(X),1);
else
    R = zeros(length(X),1);
end
for i=1:length(X)
     temp = D1(i,:);
     [x2,y2,z2]=find(temp);
     R(i)= min(z2)/2;
end

cube = (NF+1).^3;
%dist = dist./2;
sphere = 0;
for i =1:length(X)
    sphere = sphere + (4/3).*pi.*(R(i).^3);
end
density = (sphere./cube);
diameter = mean(R(:))*2;
end
