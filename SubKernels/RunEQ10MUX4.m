function [G] = RunEQ10MUX4(NF,B)
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

impares = 1:2:NF;

% Crear el vector de números pares: empieza en 2, salta de 2 en 2 hasta B
pares = 2:2:NF;

% Concatenar ambos vectores para obtener el resultado final
muxvect = [impares, pares];


for a=1:K
    for b=a:K
        for c=1:K
            for k=muxvect(1:B)
                G(:,:,k) = (mod(a*I+b*J+((k-1)*c),NF))+1;
            end
            [d] = Computedistance(G,NF);
            distance1(a,b,c) = d; %sphere diameter
        end
    end
end

[G,diameter,a1,b1,c1] = bestPattern(NF,B,distance1);
%G=G(:,:,1:2:end);
end