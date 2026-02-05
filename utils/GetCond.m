function [conditionalNumber] = GetCond(C)
NF=size(C,3);
sigma = [];
for j=1:NF
    temp = C(:,:,j);
    sigma = [sigma temp(:)];
end
sigma=sigma/max(sigma(:));
conditionalNumber = cond(sigma'*sigma);
% matrix=sigma'*sigma;
% [U,S,V] = svd(matrix);
% s=diag(S);
% val=unique(s);
% if min(val)<=0
%     disp("Warning in CN")
% end
% val(val<=0)=nan;
% conditionalNumber = max(val)/min(val);
%disp("Normalized Codition Number "+ conditionalNumber)
end