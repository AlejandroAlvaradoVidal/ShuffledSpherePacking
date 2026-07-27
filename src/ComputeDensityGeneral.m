function [density,diameter]= ComputeDensityGeneral(NF,G)
pat=ones(size(G));
for i=1:NF
    pat(:,:,i)=i.*pat(:,:,i);
end
G=G.*pat;
[r,c,z1] = find(G);
X = [r c z1]'; % generator matrix
D = pdist(X');
D1 = squareform(D);
R = zeros(numel(G),1);
for i=1:nnz(G)
     temp = D1(i,:);
     [x2,y2,z2]=find(temp);
     R(i)= min(z2)/2;
end

cube = (NF+1).^3;
%dist = dist./2;
sphere = 0;
for i =1:numel(G)
    sphere = sphere + (4/3).*pi.*(R(i).^3);
end
density = (sphere./cube);
diameter = mean(R(:))*2;
disp("Sphere Packing Density "+ density + " diameter "+ diameter);

% nb =NF.^2;
% dis = 2*(((NF+1).^3)./(4*nb*sqrt(2))).^(1/3);
% dis = dis./2;
% sphere = (NF.^2).*(4/3).*pi.*(dis.^3);
% density1 = (sphere./cube);
% disp("Optimal Sphere Packing Density "+ density1)
end