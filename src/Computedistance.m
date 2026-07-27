function [d] = Computedistance(G,NF)
[r,c,z] = find(G(1:NF,1:NF,:));
X = [r c z]'; % generator matrix
D = pdist(X');
d = min(D);
end