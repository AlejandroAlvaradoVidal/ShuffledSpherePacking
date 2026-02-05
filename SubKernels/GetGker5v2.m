function [C,G]=GetGker5v2(NF,M)
randomNumber=6;
rng(randomNumber);
for i=1:M*NF
    A(i,:) = randperm(NF);
end
C=zeros(NF,NF,NF);
c=1;
for i=1:NF:M*NF
    C=C+G2C(A(i:i+NF-1,:));
    G(:,:,c)=A(i:i+NF-1,:);
    c=c+1;
end
end