function [C] = RunEQ10MUX(NF,B)
N = B*NF; % Detector size
x = 1:NF;
y = x';
g = ones(NF,1);
q = (1:NF)';
I = kron(g',q);
J = I';
K = floor(NF/2);
distance1 = zeros(NF,NF,NF);
G = zeros(NF,NF,B);

for a=1:K
    for b=a:K
        for c=1:K
            for k=1:B
                G(:,:,k) = (mod(a*I+b*J+((k-1)*c),NF))+1;
            end
            [d] = Computedistance(G,NF);
            distance1(a,b,c) = d; %sphere diameter
        end
    end
end

[G,diameter,a1,b1,c1] = bestPattern(NF,B,distance1);
[C,T] = generateCodedAperture(G,NF,B,NF); 
end