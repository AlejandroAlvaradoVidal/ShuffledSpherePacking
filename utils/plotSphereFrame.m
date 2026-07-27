function plotSphereFrame(C, R_mode, iter, M_val, current_cond, current_density, export_path)
    % Crea y guarda un frame de la visualización de esferas para un video.
    %
    % Inputs:
    % - C: Matriz de configuración de esferas 3D.
    % - R_mode: Modo de radio (1 para irregular, 0 para regular).
    % - iter: Número de la iteración actual.
    % - M_val: Valor de Mux actual.
    % - current_cond: Valor del 'condition number' en la iteración actual.
    % - current_density: Valor de la densidad en la iteración actual.
    % - export_path: Directorio donde se guardará el frame.

    % Crear una figura invisible para no interrumpir el proceso
    fig_video = figure('Name', ['Iter ' num2str(iter)], 'Visible', 'off', 'Position', [100, 100, 800, 800]);
    
    % Utilizar la función existente para dibujar la configuración de esferas
    showSphereCA(C, 12, R_mode); 
    
    % Añadir un título informativo a la gráfica
    title(sprintf('M=%d, Iter: %d | Cond: %.2f | Density: %.4f', M_val, iter, current_cond, current_density), 'Color', 'k', 'FontSize', 14);
    
    % Guardar la figura como una imagen (frame del video)
    frame_filename = fullfile(export_path, sprintf('frame_M%d_%05d.png', M_val, iter));
    saveas(fig_video, frame_filename);
    close(fig_video); % Cerrar la figura para liberar memoria
end