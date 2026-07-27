function [C,T] = generateCodedAperture(G,NF,B,N)

N1 = ceil(N/NF);
II = ones(N1,N1);
C = zeros(NF,NF,NF);
T = zeros(N,N,NF);

for j=1:NF
    for i=1:B
        temp = G(:,:,i)==j;
        C(:,:,j) = C(:,:,j) + temp;
    end
end

transmittance = sum(C(:))/numel(C);
%disp("Coded aperture transmittance "+ transmittance);
sigma = [];

for j=1:NF
    temp = C(:,:,j);
    sigma = [sigma temp(:)];
end
sigma=sigma/(max(sigma(:)));
conditionNumber = cond(sigma'*sigma);
%disp("Codition Number "+ conditionNumber)

for i=1:NF
    temp = kron(II,C(:,:,i));
    T(:,:,i) = temp(1:N,1:N);
end
end