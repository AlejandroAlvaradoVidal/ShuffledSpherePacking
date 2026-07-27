function [C]=GetGker4(NF,M)
randomNumber=6;
%randomNumber=randi([1,100],1,1);
N=NF;
rng(randomNumber);
for i=1:NF
    A(i,:) = randperm(NF);
end
G(:,:,1) = A;

for i=2:M
    G(:,:,i)=circshift(G(:,:,i-1),1,2);
end
C=zeros(size(A));
for i=1:M
    C=C+G2C(G(:,:,i));
end
end