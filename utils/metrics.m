function [p,s,r,sam] = metrics(data,Xrec) %1 for cave dataset
[M,N,L] = size(Xrec);
p = zeros(L,1);
s = zeros(L,1);
rm = zeros(L,1);
data=mat2gray(data);
Xrec=mat2gray(Xrec);
for i=1:L
    ref = data(:,:,i);
    A = Xrec(:,:,i);
    p(i) = psnr(A,ref);
    s(i) = ssim(A,ref);
    v1 = double(A(:));
    v2 = double(ref(:));
    rm(i) = sqrt(immse(v1(:),v2(:)));
end

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

p = mean(p); % psnr
s= mean(s); % ssim
r = mean(rm); % rmse
sam = mean(sa(:)); % sam Measure spectral similarity using spectral angle mapper
disp("PSNR "+ num2str(p)+ " SSIM "+num2str(s)+" RMSE "+num2str(r)+" SAM "+num2str(sam));
disp("---------------------------------------------------------------------------------------------------------")
end