function [C]=GetGker7(NF,M)
randomNumber=6;
%randomNumber=randi([1,100],1,1);
N=NF;
rng(randomNumber);
A=rand(NF,NF,NF);
A(A<M/NF)=0;
A(A>=M/NF)=1;
C=1-A;
end