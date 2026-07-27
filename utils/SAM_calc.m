function [sam] = SAM_calc(data,Xrec) %1 for cave dataset
[M,N,L] = size(Xrec);
sa = zeros(M,N);
for n=1:N
    for m=1:M
        v1 = Xrec(m,n,:);
        v2 = data(m,n,:);
        v1 = double(v1(:));
        v2 = double(v2(:));
        sa(m,n) = real(SpectralAngleMapper(v1,v2+eps));
    end
end

sam = mean(sa(:)); 