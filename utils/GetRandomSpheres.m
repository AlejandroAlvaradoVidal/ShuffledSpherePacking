function P = GetRandomSpheres(nWant, Width, Radius)
% INPUT:
%   nWant:  Number of spheres
%   Width:  Dimension of 3d box as [1 x 3] double vector
%   Radius: Radius of spheres
% OUTPUT:
%   P:      [nWant x 3] matrix, centers

P = zeros(nWant, 3);

R2     = (2 * Radius) ^ 2;   % Squared once instead of SQRT each time
W      = Width - 2 * Radius; % Avoid interesction with borders
iLoop  = 1;                  % Security break to avoid infinite loop
nValid = 0;
while nValid < nWant && iLoop < 1e6
  newP = rand(1, 3) .* W + Radius;
  % Auto-expanding, need Matlab >= R2016b. For earlier versions:
  % Dist2 = sum(bsxfun(@minus, P(1:nValid, :), newP) .^ 2, 2);
  Dist2 = sum((P(1:nValid, :) - newP) .^ 2, 2);
  if all(Dist2 > R2)
    % Success: The new point does not touch existing sheres:
    nValid       = nValid + 1;  % Append this point
    P(nValid, :) = newP;
  end
  iLoop = iLoop + 1;
end
% Stop if too few values have been found:
if nValid < nWant
  error('Cannot find wanted number of points in %d iterations.', iLoop)
end
end