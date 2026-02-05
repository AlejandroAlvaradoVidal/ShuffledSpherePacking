function [val] = findhip(x,M)
load('hip.mat')
y = (x/M).*ones(size(hip));
[~, idx] = min(abs(hip - y));

% Obtén el valor aproximado
val = hip(idx);

end