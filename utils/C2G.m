function [G] = C2G(C)
[N,M,B]=size(C);

G=zeros(N,M);
for k=1:B
    for i=1:N
        for j=1:M
            if C(i,j,k)==1
                G(i,j)=k;
            end
        end
    end
end
end