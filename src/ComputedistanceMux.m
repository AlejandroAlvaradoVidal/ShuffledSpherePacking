function [d] = ComputedistanceMux(G,NF)
%G=G(1:NF,1:NF,:);
Gt=[];
if size(G,3)~=1
    for i=1:size(G,3)
    Gt=vertcat(Gt,G(:,:,i));
    end
else
    Gt=G;
end
[r,c,z] = find(Gt);
X = [r c z]'; % generator matrix
D = pdist(X');
d = min(D);
end