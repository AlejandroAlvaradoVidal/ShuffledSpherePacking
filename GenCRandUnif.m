addpath(genpath('./src'));
addpath(genpath('./utils'));
addpath(genpath('./Method'));
addpath(genpath('./SubKernels'));
addpath(genpath('./PlotSpheres'));
clc; clear; close all;

%% --- Parámetros de Configuración ---
NF = 31;                % Bandas (Profundidad Z)
N = 31;                 % Tamaño del Sensor
R = 1;                  % 0 para Regular, 1 para Irregular
max_iterations = 2000;  % Máximo de iteraciones
rel_tol = 1e-5;         % Mejora mínima relativa
patience = 80;          % Iteraciones permitidas sin mejora
burn_in = 10;          % Mínimo de iteraciones antes de evaluar salida

%% --- Parámetros de Visualización ---
% 0: No graficar nada (más rápido)
% 1: Exportar gráfico a carpeta sin mostrar ventana (Headless)
% 2: Gráfico tradicional (Mostrar ventana)
plot_mode = 1; 
export_path = './graphsU'; % Directorio para guardar si plot_mode = 1
sf=1; % 1 for save C

for M = 1:16
    % --- Inicialización de la Matriz ---
    %[C,G]=GetGker5v2(NF,M);
    [H, C] = GenerateSensingMatrix(16, 8);
    C(C > 1) = 1;
    
    Cond = GetCond(C);
    best_cond = Cond;
    best_C = C; % Si M=1, best_C será la generada aquí arriba
    
    % Solo optimizar y graficar si M > 1
    % Históricos para gráficas
    CNr = nan(1, max_iterations); 
    sumden = nan(1, max_iterations);
    
    % Control de salida
    no_improvement_count = 0;
    
    fprintf('Iniciando Optimización M=%d | Cond Inicial: %.2f\n', M, Cond);
    fprintf('------------------------------------------------------------\n');

    %% --- Randomizer Section (XZ, YZ, XY) ---
    for k = 1:max_iterations
        % Identificar posiciones críticas basándose en distancias
        [list, diameters] = uniqueDistanceIrregularSP(C);
        
        % Extraer coordenadas actuales de los elementos activos (1s)
        X = [];
        for i_layer = 1:NF
            tp = C(:,:,i_layer) * i_layer;
            [x_idx, y_idx, z_idx] = find(tp);
            X = [X; [x_idx, y_idx, z_idx]];
        end
        
        % Seleccionamos un punto aleatorio dentro de los diámetros mínimos
        positions = find(diameters == list(1));
        val = randi([1, length(positions)], 1, 1);
        [xyzt] = X(positions(val), :);
        i = xyzt(1); ii = xyzt(2); j = xyzt(3);
        
        % Copia temporal para evaluar el movimiento
        C_next = C;
        
        % --- Selección de Plano de Permutación ---
        selector = rand();
        
        if selector < 0.5
            % --- INTERCAMBIO PLANO XZ ---
            otherLayer = randi([1, NF]);
            temp = C(i, :, j); 
            C(i, :, j) = C(i, :, otherLayer);
            C(i, :, otherLayer) = temp;
        else
            % --- INTERCAMBIO PLANO YZ ---
            otherLayer = randi([1, NF]);
            temp = C(:, ii, j);
            C(:, ii, j) = C(:, ii, otherLayer);
            C(:, ii, otherLayer) = temp;
        end
        
        C_next(C_next > 1) = 1;
        
        % --- Evaluación de la Mejora ---
        new_cond = GetCond(C_next);
        if R == 1
            [new_density, ~] = ComputeDensityIrregularSP(C_next);
        else
            [new_density, ~] = ComputeDensityRegularSP(C_next);
        end

        if new_cond < best_cond
            improvement = (best_cond - new_cond) / best_cond;
            if improvement > rel_tol
                best_cond = new_cond;
                best_C = C_next;
                C = C_next; 
                no_improvement_count = 0;
            else
                no_improvement_count = no_improvement_count + 1;
            end
        else
            no_improvement_count = no_improvement_count + 1;
        end

        CNr(k) = new_cond;
        sumden(k) = new_density;

        if mod(k, 20) == 0
            fprintf('Iter: %d | Cond: %.2f | Best: %.2f | Status: %d/%d\n', ...
                k, new_cond, best_cond, no_improvement_count, patience);
        end

        if k > burn_in && no_improvement_count >= patience
            fprintf('>> Convergencia alcanzada en iteración %d.\n', k);
            break;
        end
    end

    %% --- Generación de Gráficos ---
    if plot_mode > 0
        % Configuración de visibilidad según el modo
        fig_vis = 'on';
        if plot_mode == 1, fig_vis = 'off'; end
        
        % Crear la figura con la visibilidad deseada
        fig = figure('Name', ['Optimización Mux ' num2str(M)], ...
                        'Color', 'w', ...
                        'Position', [100, 100, 800, 600], ...
                        'Visible', fig_vis);
        
        % --- Subplot 1: Condition Number ---
        subplot(2,1,1);
        plot(1:k, CNr(1:k), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
        hold on; 
        yline(best_cond, '--r', ['Min: ' num2str(best_cond, '%.2f')], 'LabelVerticalAlignment', 'bottom');
        
        set(gca, 'YScale', 'log', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k');
        box on; title(['Evolución del Condition Number (M=' num2str(M) ')'], 'Color', 'k');
        ylabel('Cond(C)'); xlabel('Iteraciones');
        
        % --- Subplot 2: Density ---
        subplot(2,1,2);
        plot(1:k, sumden(1:k), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
        
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k');
        box on; title('Evolución de la Densidad', 'Color', 'k');
        ylabel('Densidad'); xlabel('Iteraciones');
        
        % Lógica de exportación si el modo es 1
        if plot_mode == 1
            if ~exist(export_path, 'dir'), mkdir(export_path); end
            filename = fullfile(export_path, ['Opt_Mux_', num2str(M), '.png']);
            saveas(fig, filename);
            close(fig); % Cerramos la figura invisible para liberar RAM
            fprintf('>> Gráfico M=%d exportado en: %s\n', M, filename);
        end
    end

    %% --- Procesamiento Final de la Matriz ---
    % best_C contiene la matriz inicial si M=1, o la optimizada si M>1
    N1 = ceil(N/NF);
    II = ones(N1,N1);
    T = zeros(N, N, NF);
    for i_f = 1:NF
        temp_kron = kron(II, best_C(:,:,i_f));
        T(:,:,i_f) = temp_kron(1:N, 1:N);
    end
    C_final = T;
    C = C_final;

    fprintf('Proceso M=%d finalizado. Condición final: %.4f\n\n', M, best_cond);
    
    if sf==1
        folderName = "./Patrones/RUnif";
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        % Nota: Asegúrate que la variable 'name' esté definida o usa 'folderName'
        save(fullfile(folderName, ['Bands_', num2str(NF), '_Mux_', num2str(M), '.mat']), "C");
    end
end