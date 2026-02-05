function [G]=GetGker(NF,M)
randomNumber=6;
N=NF;
rng(randomNumber);
for i=1:NF
    A(i,:) = randperm(NF);
end
r = ceil(N/NF);
B = ones(r,r);
G = kron(B,A);
G(:,:,1) = G(1:N,1:N);
temp=G;
for i=2:M
    temp=temp(randperm(end),randperm(end));
    G(:,:,i)=temp;
end
end