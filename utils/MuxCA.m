% Created by: Ph.D Nelson Eduardo Diaz Diaz
% Post-doctorado Pontícia Universidad Católica de Valparaíso (PUCV)
% Date 2 February 2022

% Comparison of minimum distance
% Solution using the Discrete Sphere Packing based on 3D N^2 Queens
% Approach

% Find Optimal values of a and b

function [a,b,ma,G]=MuxCA(N,NF,B)
M = round(NF/2);
x = 1:NF;
y = x';
g = ones(NF,1);
q = (1:NF)';
I = kron(g',q);
J = I';
K = floor(NF/2);
distance = zeros(NF,NF,NF);
G = (zeros(NF,NF,B));

for i=1:M
    for j=i:M
        for c=1:M
            for k=1:B
                G(:,:,k) = (mod(I.*i + J.*j+((k-1)*c),NF))+1;
            end
            [d] = ComputeDensityGeneralMux(NF,G,B);
        end
        %disp(i +" out of "+ num2str(M));
    end
end
ma = max(distance(:));
[G,diameter,a1,b1,c1] = bestPattern(NF,B,d);
% [b,a,~]= find(ma==distance);
% G = (mod(I.*a(1) + J.*b(1)+((k-1)*c(1)),NF))+1;
% A = G;
% p = ceil(N/NF);
% B = ones(p,p);
% G1 = kron(B,A);
% G1 = G1(1:N,1:N);
%save("results/dist_best"+num2str(N),'distance','a','b')
end