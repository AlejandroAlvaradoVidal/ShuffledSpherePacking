function [C] = C2G(G)
[N,M]=size(G);
B=max(G(:));
C=zeros(N,M,B);
for k=1:B
    for i=1:N
        for j=1:M
            if G(i,j)==k
                C(i,j,k)=1;
            end
        end
    end
end
end