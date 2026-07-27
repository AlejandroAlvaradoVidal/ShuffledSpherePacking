function blueNoise3D = BlueNoise3D(N, M)
    % gridSize: The size of the 3D grid [X, Y, Z]
    % numPoints: Number of points to generate
    
    % Step 1: Initialize a grid of random noise
    numPoints=M/N*N^3;
    gridSize = [N,N, N];
    randomNoise = rand(gridSize);
    
    % Step 2: Apply a high-pass filter in the frequency domain
    fftNoise = fftn(randomNoise);
    [X, Y, Z] = ndgrid(-floor(gridSize(1)/2):floor(gridSize(1)/2)-1, ...
                       -floor(gridSize(2)/2):floor(gridSize(2)/2)-1, ...
                       -floor(gridSize(3)/2):floor(gridSize(3)/2)-1);
    
    % Create a radial frequency filter
    freqRadius = sqrt(X.^2 + Y.^2 + Z.^2);
    highPassFilter = freqRadius > min(gridSize) * 0.1; % Adjust cutoff for noise properties
    
    % Apply the filter
    filteredNoise = ifftn(fftn(randomNoise) .* ifftshift(highPassFilter));
    
    % Step 3: Normalize the result to get blue noise
    blueNoise3D = real(filteredNoise);
    blueNoise3D = (blueNoise3D - min(blueNoise3D(:))) / (max(blueNoise3D(:)) - min(blueNoise3D(:)));
    
    % Step 4: Threshold to get binary blue noise halftone
    [~, indices] = sort(blueNoise3D(:), 'descend');
    halftoneMask = zeros(size(blueNoise3D));
    halftoneMask(indices(1:numPoints)) = 1;
    
    % Optional: Visualize the result
    [x, y, z] = ind2sub(size(halftoneMask), find(halftoneMask));
    scatter3(x, y, z, 'filled');
    axis equal;
    title('3D Blue Noise Halftoning');
    
    % Output
    blueNoise3D = halftoneMask;
end

