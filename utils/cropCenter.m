function croppedImage = cropCenter(image, n)
    % Corta el área central de tamaño nxn de una imagen.
    %
    % Inputs:
    % - image: Imagen de entrada (puede ser en escala de grises o RGB).
    % - n: Tamaño del recorte (n x n píxeles).
    %
    % Output:
    % - croppedImage: Imagen recortada.

    % Verificar si el tamaño solicitado es válido
    [height, width, ~] = size(image); % Dimensiones de la imagen
    if n > height || n > width
        error('El tamaño de recorte n es mayor que las dimensiones de la imagen.');
    end

    % Calcular las coordenadas del área central
    startRow = round((height - n) / 2) + 1;
    startCol = round((width - n) / 2) + 1;
    endRow = startRow + n - 1;
    endCol = startCol + n - 1;

    % Recortar el área central
    croppedImage = image(startRow:endRow, startCol:endCol, :);
end
