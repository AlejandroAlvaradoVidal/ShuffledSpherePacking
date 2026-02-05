function entropy = shannonEntropy3D(points, num_bins)
    % Calcula la entropía de Shannon para una lista de puntos 3D
    % points: una matriz Nx3 donde cada fila representa un punto [x, y, z]
    % num_bins: número de divisiones en cada dimensión para el histograma 3D

    if nargin < 2
        num_bins = 10; % Valor predeterminado si no se especifica
    end

    % Crear el histograma 3D usando hist3 en las primeras dos dimensiones
    edges = {linspace(min(points(:,1)), max(points(:,1)), num_bins + 1), ...
             linspace(min(points(:,2)), max(points(:,2)), num_bins + 1), ...
             linspace(min(points(:,3)), max(points(:,3)), num_bins + 1)};
    
    % Inicializar la matriz de conteos
    counts3D = zeros(num_bins, num_bins, num_bins);

    % Calcular los conteos en 3D manualmente
    for i = 1:num_bins
        for j = 1:num_bins
            for k = 1:num_bins
                % Verificar cuántos puntos caen dentro de cada celda
                in_bin = points(:,1) >= edges{1}(i) & points(:,1) < edges{1}(i+1) & ...
                         points(:,2) >= edges{2}(j) & points(:,2) < edges{2}(j+1) & ...
                         points(:,3) >= edges{3}(k) & points(:,3) < edges{3}(k+1);
                counts3D(i,j,k) = sum(in_bin);
            end
        end
    end

    % Normalizar las cuentas para obtener probabilidades
    prob = counts3D / sum(counts3D(:));

    % Remover probabilidades cero para evitar log2(0)
    prob(prob == 0) = [];

    % Calcular la entropía de Shannon
    entropy = -sum(prob .* log2(prob));
end
