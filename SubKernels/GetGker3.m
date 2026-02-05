function [C]=GetGker3(NF,M)
randomNumber=randi([1,NF],1,1);
rng(randomNumber);
N = NF; % Detector .
T = M/NF; %Transmittance 
m = repelem([1 0], [round(T*NF^3) round((1-T)*NF^3)]);
C = reshape(m(randperm(numel(m))), NF,NF,NF);
end