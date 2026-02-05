function [G]=GetGkerM(NF,M)
%randomNumber=6;
randomNumber=randi([1,100],1,M);
N=NF;

for j=1:M
    rng(randomNumber(j));
    for i=1:NF
        A(i,:) = randperm(NF);
    end
    G(:,:,j) = A;
end
end