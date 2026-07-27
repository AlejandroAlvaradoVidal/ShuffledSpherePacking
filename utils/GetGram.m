function [u] = GetGram(C)
NF=size(C,3);
sigma = [];
for j=1:NF
    temp = C(:,:,j);
    sigma = [sigma temp(:)];
end
sigma=sigma/max(sigma(:));
u=sigma'*sigma;
end