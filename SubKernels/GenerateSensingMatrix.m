function [H, masks] = GenerateSensingMatrix(NF, n)
    % NF: Dimensión espacial/spectral (tamaño de la máscara NF x NF)
    % n:  Número de veces que cada píxel está "activo" a través de las máscaras
    
    % 1. Inicializar la matriz H de dimensión NF x NF^2
    % Cada fila de H representará una máscara de apertura (shot)
    H = zeros(NF, NF^2);
    
    % 2. Llenar columna por columna para garantizar que la suma sea n
    for j = 1:NF^2
        % Crear un vector con n unos y el resto ceros
        columna = [ones(n, 1); zeros(NF - n, 1)];
        
        % Desordenar aleatoriamente la columna
        H(:, j) = columna(randperm(NF));
    end
    
    % 3. Extraer las máscaras individuales (códigos de apertura)
    masks = zeros(NF, NF, NF);
    for i = 1:NF
        % Reshape de la fila i para obtener una matriz NF x NF
        masks(:, :, i) = reshape(H(i, :), [NF, NF]);
    end
    
    % 4. Verificación
    suma_total = sum(masks, 3);
    fprintf('Verificación: El valor único en la suma de máscaras es: %d\n', unique(suma_total));
end