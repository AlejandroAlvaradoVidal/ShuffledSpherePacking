addpath(genpath('./src'));
addpath(genpath('./utils'));
addpath(genpath('./Method'));
addpath(genpath('./SubKernels'));
addpath(genpath('./PlotSpheres'));
clc; clear; close all;

%% --- Parámetros de Configuración ---
NF = 31;                % Bandas (Profundidad Z)
N = NF;                 % Tamaño del Sensor
R = 1;                  % 0 para Regular, 1 para Irregular
max_iterations = 2000;  % Máximo de iteraciones
rel_tol = 1e-5;         % Mejora mínima relativa
patience = 80;          % Iteraciones permitidas sin mejora
burn_in = 10;           % Mínimo de iteraciones antes de evaluar salida

%% --- Parámetros de Visualización ---
% 0: No graficar nada (más rápido)
% 1: Exportar gráfico a carpeta sin mostrar ventana (Headless)
% 2: Gráfico tradicional (Mostrar ventana)
plot_mode = 1; 
export_path = './graphsv3'; % Directorio para guardar si plot_mode = 1
sf=1; % 1 for save C

for M = 1:16
    % --- Inicialización de la Matriz ---
    % Se genera una matriz C con 'unos' en posiciones totalmente aleatorias.
    % El número total de 'unos' es M * NF, que es el esperado.
    C = zeros(NF, NF, NF);
    total_elements = NF * NF * NF;
    elements_to_set = M * NF;
    indices = randperm(total_elements, elements_to_set);
    C(indices) = 1;
    C(C > 1) = 1;
    
    Cond = GetCond(C);
    best_cond = Cond;
    best_C = C;
    
    if M > 1
        % Históricos para gráficas
        CNr = nan(1, max_iterations); 
        sumden = nan(1, max_iterations);
        
        no_improvement_count = 0;
        
        fprintf('Iniciando Optimización M=%d | Cond Inicial: %.2f\n', M, Cond);
        fprintf('------------------------------------------------------------\n');

        flag25 = 0; % Inicializar flag
        if M/NF < 0.25
            flag25 = 1; 
            fselector = rand();
        end

        %% --- Bucle de Optimización ---
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
            
            if isempty(list) || isempty(X)
                fprintf('>> No se encontraron puntos críticos. Saltando iteración.\n');
                continue;
            end

            % Seleccionamos un punto aleatorio dentro de los diámetros mínimos
            positions = find(diameters == list(1));
            val = randi([1, length(positions)], 1, 1);
            xyzt = X(positions(val), :);
            i = xyzt(1); ii = xyzt(2); j = xyzt(3);
            
            C_next = C;
                        
            % --- Selección de Plano de Permutación (Lógica de v2) ---
            use_xz_plane = false;
            if flag25 == 1
                use_xz_plane = (fselector < 0.5);
            else
                use_xz_plane = (rand() < 0.5);
            end

            otherLayer = randi([1, NF]);
            if use_xz_plane
                % --- INTERCAMBIO PLANO XZ ---
                temp = C_next(i, :, j); 
                C_next(i, :, j) = C_next(i, :, otherLayer);
                C_next(i, :, otherLayer) = temp;
            else
                % --- INTERCAMBIO PLANO YZ ---
                temp = C_next(:, ii, j);
                C_next(:, ii, j) = C_next(:, ii, otherLayer);
                C_next(:, ii, otherLayer) = temp;
            end
           
            C_next(C_next > 1) = 1;

            % --- PROCESO DE UNIFORMACIÓN ---
            C_next = uniformizeMatrix(C_next, M);
            
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
            fig_vis = 'on';
            if plot_mode == 1, fig_vis = 'off'; end
            
            fig = figure('Name', ['Optimización Mux ' num2str(M)], 'Color', 'w', 'Position', [100, 100, 800, 600], 'Visible', fig_vis);
            
            subplot(2,1,1);
            plot(1:k, CNr(1:k), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
            hold on; 
            yline(best_cond, '--r', ['Min: ' num2str(best_cond, '%.2f')], 'LabelVerticalAlignment', 'bottom');
            set(gca, 'YScale', 'log', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 1);
            grid on;
            box on; title(['Evolución del Condition Number (M=' num2str(M) ')'], 'Color', 'k');
            ylabel('Cond(C)'); xlabel('Iteraciones');
            
            subplot(2,1,2);
            plot(1:k, sumden(1:k), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
            set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 1);
            grid on;
            box on; title('Evolución de la Densidad', 'Color', 'k');
            ylabel('Densidad'); xlabel('Iteraciones');
            
            if plot_mode == 1
                if ~exist(export_path, 'dir'), mkdir(export_path); end
                filename = fullfile(export_path, ['Opt_Mux_', num2str(M), '.png']);
                saveas(fig, filename);
                close(fig);
                fprintf('>> Gráfico M=%d exportado en: %s\n', M, filename);
            end
        end
    else
        fprintf('M=1 detected: Skip Method.\n');
    end

    %% --- Procesamiento Final de la Matriz ---
    N1 = ceil(N/NF);
    II = ones(N1,N1);
    T = zeros(N, N, NF);
    for i_f = 1:NF
        temp_kron = kron(II, best_C(:,:,i_f));
        T(:,:,i_f) = temp_kron(1:N, 1:N);
    end
    C_final = T;

    fprintf('Proceso M=%d finalizado. Condición final: %.4f\n\n', M, best_cond);
    
    if sf==1
        folderName = "./Patrones/SSPv3";
        if ~exist(folderName, 'dir'), mkdir(folderName); end
        save(fullfile(folderName, ['Bands_', num2str(NF), '_Mux_', num2str(M), '.mat']), "C_final");
    end
end

function C_out = uniformizeMatrix(C_in, M)
    % Ajusta la matriz C para que la suma en cada dimensión sea M.
    C_out = C_in;
    [dim1, dim2, dim3] = size(C_out);

    for dim = 1:3 % Itera sobre las 3 dimensiones
        for i = 1:size(C_out, dim)
            % Extrae una "rebanada" 2D
            if dim == 1, slice_indices = {i, ':', ':'};
            elseif dim == 2, slice_indices = {':', i, ':'};
            else, slice_indices = {':', ':', i};
            end
            
            slice = C_out(slice_indices{:});
            current_sum = sum(slice, 'all');
            diff = current_sum - M;

            if diff > 0 % Hay 'unos' de más
                ones_idx_flat = find(slice);
                % Asegurarse de no intentar eliminar más 'unos' de los que existen
                num_to_move = min(diff, length(ones_idx_flat));
                ones_to_move = ones_idx_flat(randperm(length(ones_idx_flat), num_to_move));
                slice(ones_to_move) = 0;
                C_out(slice_indices{:}) = slice;

            elseif diff < 0 % Faltan 'unos'
                zeros_idx_flat = find(~slice);
                % Asegurarse de no intentar añadir más 'unos' de los que hay espacio
                num_to_fill = min(-diff, length(zeros_idx_flat));
                zeros_to_fill = zeros_idx_flat(randperm(length(zeros_idx_flat), num_to_fill));
                slice(zeros_to_fill) = 1;
                C_out(slice_indices{:}) = slice;
            end
        end
    end
end