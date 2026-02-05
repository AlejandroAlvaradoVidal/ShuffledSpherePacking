% Created by: Ph.D Nelson Eduardo Diaz Diaz
% Post-doctorado Pontícia Universidad Católica de Valparaíso (PUCV)
% Date 2 February 2022

% Comparison of minimum distance
% Solution using the Discrete Sphere Packing based on 3D N^2 Queens
% Approach

% Find Optimal values of a and b

function [a,b,ma,G1]=DDDRSNNP3(N,NF)
M = round(NF/2);
distance = zeros(M,M);
x = ones(1,NF)';
y = (1:NF)';
I = kron(x',y);
J = kron(x,y');

for i=1:M
    for j=i:M
        G = mod(I.*i + J.*j,NF)+1;
        t = length(unique(G(1:NF,1:NF)));
        if(t == NF)
            G1 = G(1:NF,1:NF);
            [r,c,z] = find(G1);
            B = [r c z];
            D = pdist(B);
            [distance(i,j)]=min(D(:));
        else
        end
    end
    %disp(i +" out of "+ num2str(M));
end
ma = max(distance(:));
[b,a,~]= find(ma==distance);
G = mod(I.*a(1) + J.*b(1),NF)+1;
A = G;
p = ceil(N/NF);
B = ones(p,p);
G1 = kron(B,A);
G1 = G1(1:N,1:N);
%save("results/dist_best"+num2str(N),'distance','a','b')
end