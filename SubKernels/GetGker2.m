function [C]=GetGker2(NF,M)
randomNumber=6;%randi([1,NF],1);
for i=1:NF
    A(i,:) = randperm(NF);
end
G(:,:,1) = A;
% [a,b,ma,G1]=DDDRSNNP3(NF,M)
 G(:,:,1) = A;
 temp=G;
for i=2:M
    temp=G(:,:,i-1);
    temp=temp(randperm(end),randperm(end));
    G(:,:,i)=temp;
end
C=zeros(NF,NF,NF);
for i=1:M
    C=C+G2C(G(:,:,i));
end

end