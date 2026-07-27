function [conditionalNumber] = GetCond(C)
NF=size(C,3);
is_gpu = isa(C, 'gpuArray');
if is_gpu
    sigma = gpuArray([]);
else
    sigma = [];
end

for j=1:NF
    temp = C(:,:,j);
    sigma = [sigma temp(:)];
end
sigma=sigma/max(sigma(:));
conditionalNumber = cond(sigma'*sigma);
end