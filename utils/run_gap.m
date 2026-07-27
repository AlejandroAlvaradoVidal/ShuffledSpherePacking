function [vgaptv] = run_gap(orig,mask)
addpath(genpath('./algorithms')); % algorithms
addpath(genpath('./packages')); % packages
addpath(genpath('./utils_gap')); % utilities

meas=sum(orig.*mask,3);

para.nframe =   1; % number of coded frames in this test
para.MAXB   = 255;

[nrow,ncol,nmask] = size(mask);
nframe = para.nframe; % number of coded frames in this test
MAXB = para.MAXB;

% [1.2] parameter setting for GAP-TV and GAP-WNNM
para.Mfunc  = @(z) A_xy(z,mask);
para.Mtfunc = @(z) At_xy_nonorm(z,mask);

para.Phisum = sum(mask.^2,3);
para.Phisum(para.Phisum==0) = 1;
% common parameters
para.lambda   =     1; % correction coefficiency
para.acc      =     1; % enable GAP-acceleration
para.flag_iqa = false; % disable image quality assessments in iterations

%% [2.1] GAP-TV, ICIP'16
para.lambda   =    1; % correction coefficiency
para.maxiter  =  300; % maximum iteration
para.acc      =    1; % enable acceleration
para.denoiser = 'tv'; % TV denoising
  para.tvweight = 0.07*255/MAXB; % weight for TV denoising
  para.tviter   = 5; % number of iteration for TV denoising
  
[vgaptv,psnr_gaptv,ssim_gaptv,tgaptv] = ...
    gapdenoise_cacti(mask,meas,orig,[],para);

end