function [density,diameter]= ComputeDensityGeneralMux(NF,C,B)
pat=ones(size(C));
for i=1:NF
    pat(:,:,i)=i.*pat(:,:,i);
end
x=C.*pat;
[r,c,z1] = find(x);
X = [r c z1]'; % generator matrix
D = pdist(X');
D1 = squareform(D);
R = zeros(nnz(C),1);
for i=1:nnz(C)
     temp = D1(i,:);
     [x2,y2,z2]=find(temp);
     R(i)= min(z2)/2;
end

cube = (NF+1).^3;
%dist = dist./2;
sphere = 0;
for i =1:nnz(C)
    sphere = sphere + (4/3).*pi.*(R(i).^3);
end
density = (sphere./cube);
diameter = mean(R(:))*2;
end