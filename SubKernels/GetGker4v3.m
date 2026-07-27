function [C]=GetGker4v3(NF,M)
randomNumber=6;
%randomNumber=randi([1,100],1,1);
N=NF;
rng(randomNumber);
for i=1:M*NF
    A(:,:,i) = randperm(NF,NF);
end

for i=1:M
    G(:,:,i) = A(:,:,i);
end

C=zeros(NF,NF,NF);
for i=1:M
    C=C+G2C(G(:,:,i));
end
end