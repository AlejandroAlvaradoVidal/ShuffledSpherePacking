 function [G,diameter,a1,b1,c1] = bestPattern(NF,B,MDistance)

g = ones(NF,1);
q = (1:NF)';
I = kron(g',q);
J = I';

[maxval,ind] = max(reshape(MDistance(:),NF*NF*NF,[]));
[a1,b1,c1]= ind2sub([NF NF NF],ind);
diameter = maxval;

for k=1:B
    G(:,:,k) = (mod(a1(1)*I+b1(1)*J+((k-1)*c1(1)),NF))+1;
end
end