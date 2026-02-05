% Created by: Ph.D Nelson Eduardo Diaz Diaz
% Post-doctorado Pontícia Universidad Católica de Valparaíso (PUCV)
% Date 2 February 2022

% Amidror, Isaac. “Scattered data interpolation methods for electronic imaging systems: a survey.” 
% Journal of Electronic Imaging. Vol. 11, No. 2, April 2002, pp. 157–176.

function [Xrec] = interpolation(J)
[M,N,L]= size(J);
Xrec = zeros(M,N,L);
md = 5;
for j=1:L
    X = J(:,:,j);
    [x,y,v] = find(X);
    [xq,yq] = find(X==0);
    %x1 = interp2(x,y,v,xq,yq);
    F = scatteredInterpolant(x,y,v,'natural');
    X(X==0) = F(xq,yq);
    Xrec(:,:,j) = X; 
    %Xrec(:,:,j) = medfilt2(Xrec(:,:,j),[md md]);
    Xrec(:,:,j) = imsharpen(Xrec(:,:,j));
end
end