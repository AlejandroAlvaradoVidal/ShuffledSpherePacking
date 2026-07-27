addpath(genpath('./src'));
addpath(genpath('./utils'));
addpath(genpath('./Method'));
addpath(genpath('./SubKernels'));
addpath(genpath('./PlotSpheres'));
clc; clear; close all;

%% --- Parámetros de Configuración ---
NF = 31;                % Bands (Z-depth)
N = NF;                 % Sensor Size
R = 1;                  % 0 for Regular, 1 for Irregular
max_iterations = 2000;  % Maximum iterations
rel_tol = 1e-5;         % Minimum relative improvement
patience = 2000;          % Iterations allowed without improvement
optimization_metric = 'density'; % 'cond' for condition number, 'density' for maximizing density
burn_in = 10;           % Minimum iterations before evaluating exit

%% --- Parámetros de Visualización ---
% 0: Do not plot anything (faster)
% 1: Export graph to folder without showing window (Headless)
% 2: Traditional plot (Show window)
plot_mode = 1; 
show_sphere_plots = true; % Bypass for sphere visualization
export_path = './graphsv4'; % Directory to save if plot_mode = 1
sf=1; % 1 for save C

% --- Configuración de Archivo de Log ---
log_filename = 'summary_log.txt';
fileID = fopen(log_filename, 'w'); % Open in write mode to clear on each run
fprintf(fileID, 'Execution Summary - %s\n\n', datetime('now'));

mux_values_to_run = 1:16;

% --- Arrays para almacenar métricas para el gráfico de resumen ---
summary_initial_cond = nan(1, length(mux_values_to_run));
summary_optimized_cond = nan(1, length(mux_values_to_run));
summary_uniform_cond = nan(1, length(mux_values_to_run));
summary_initial_density = nan(1, length(mux_values_to_run));
summary_optimized_density = nan(1, length(mux_values_to_run));
summary_uniform_density = nan(1, length(mux_values_to_run));
summary_initial_rmse = nan(1, length(mux_values_to_run));
summary_optimized_rmse = nan(1, length(mux_values_to_run));
summary_uniform_rmse = nan(1, length(mux_values_to_run));

for m_idx = 2:length(mux_values_to_run)
    M = mux_values_to_run(m_idx);
    % --- Inicialización de la Matriz ---
    [G] = RunEQ10MUX5(NF, M);
    [C, ~] = generateCodedAperture(G, NF, M, NF); 
    C(C > 1) = 1;
    C_initial = C; % Guardar matriz inicial
    
    % Calcular métricas iniciales
    initial_cond = GetCond(C);
    % Corrección para M=1: El número de condición de una matriz de permutación ideal es 1.
    if M == 1
        initial_cond = 1;
    end

    if R == 1
        [initial_density, ~] = ComputeDensityIrregularSP(C);
    else
        [initial_density, ~] = ComputeDensityRegularSP(C);
    end

    % Calculate initial uniformity RMSE
    sum_x_init = reshape(sum(C_initial, [2, 3]), [N, 1]);
    sum_y_init = reshape(sum(C_initial, [1, 3]), [N, 1]);
    sum_z_init = reshape(sum(C_initial, [1, 2]), [N, 1]);
    errors_init = [(sum_x_init - M); (sum_y_init - M); (sum_z_init - M)];
    initial_rmse = sqrt(mean(errors_init.^2));

    best_cond = initial_cond;
    best_density = initial_density;
    density_after = initial_density; % Initialize for M=1 case
    best_C = C; % If M=1, best_C will be the one generated above
    final_rmse = initial_rmse; % Initialize for M=1 case
    
    % Solo optimizar y graficar si M > 1
    if M > 1
        % History for plots
        CNr = nan(1, max_iterations); 
        sumden = nan(1, max_iterations);
        
        % Exit control
        no_improvement_count = 0;
        
        fprintf('Starting Optimization M=%d | Initial Cond: %.2f | Initial Density: %.4f\n', M, initial_cond, initial_density);
        fprintf('------------------------------------------------------------\n');

        flag25 = 0; % Initialize flag
        if M/NF < 0.25
            flag25 = 1; 
            fselector = rand();
        end
        %% --- Optimization Loop (v2 Logic) ---
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
                fprintf('>> No critical points found. Skipping iteration.\n');
                continue;
            end

            % Seleccionamos un punto aleatorio dentro de los diámetros mínimos
            positions = find(diameters == list(1));
            val = randi([1, length(positions)], 1, 1);
            xyzt = X(positions(val), :);
            i = xyzt(1); ii = xyzt(2); j = xyzt(3);
            
            C_next = C;
                        
            % --- Permutation Plane Selection (Refactored) ---
            use_xz_plane = false;
            if flag25 == 1
                use_xz_plane = (fselector < 0.5);
            else
                use_xz_plane = (rand() < 0.5);
            end

            otherLayer = randi([1, NF]);
            if use_xz_plane
                % --- XZ PLANE SWAP ---
                temp = C_next(i, :, j); 
                C_next(i, :, j) = C_next(i, :, otherLayer);
                C_next(i, :, otherLayer) = temp;
            else
                % --- YZ PLANE SWAP ---
                temp = C_next(:, ii, j);
                C_next(:, ii, j) = C_next(:, ii, otherLayer);
                C_next(:, ii, otherLayer) = temp;
            end
           
            C_next(C_next > 1) = 1;
           
            % --- Improvement Evaluation ---
            new_cond = GetCond(C_next);
            if R == 1
                [new_density, ~] = ComputeDensityIrregularSP(C_next);
            else
                [new_density, ~] = ComputeDensityRegularSP(C_next);
            end

            % --- Acceptance Logic (configurable) ---
            improvement_found = false;
            if strcmpi(optimization_metric, 'cond')
                % Objective: minimize condition number
                if new_cond < best_cond
                    is_significant = isinf(best_cond) || ((best_cond - new_cond) / best_cond > rel_tol);
                    if is_significant
                        best_cond = new_cond;
                        best_density = new_density; % Update secondary metric
                        improvement_found = true;
                    end
                end
            elseif strcmpi(optimization_metric, 'density')
                % Objective: maximize density
                if new_density > best_density
                    % Avoid division by zero if best density is 0
                    if best_density > 0
                        improvement_ratio = (new_density - best_density) / best_density;
                    else
                        improvement_ratio = inf; % Any density > 0 is an infinite improvement over 0
                    end
                    
                    if improvement_ratio > rel_tol
                        best_density = new_density;
                        best_cond = new_cond; % Actualizar métrica secundaria
                        improvement_found = true;
                    end
                end
            end

            if improvement_found
                best_C = C_next;
                C = C_next; 
                no_improvement_count = 0;
            else
                no_improvement_count = no_improvement_count + 1;
            end

            CNr(k) = new_cond;
            sumden(k) = new_density;

            if mod(k, 20) == 0
                if strcmpi(optimization_metric, 'cond')
                    fprintf('Iter: %d | Cond: %.2f | Best Cond: %.2f | Status: %d/%d\n', k, new_cond, best_cond, no_improvement_count, patience);
                else % 'density'
                    fprintf('Iter: %d | Dens: %.4f | Best Dens: %.4f | Status: %d/%d\n', k, new_density, best_density, no_improvement_count, patience);
                end
            end

            if k > burn_in && no_improvement_count >= patience
                fprintf('>> Convergence reached at iteration %d.\n', k);
                break;
            end
        end

        %% --- Plot Generation ---
        if plot_mode > 0
            fig_vis = 'on';
            if plot_mode == 1, fig_vis = 'off'; end
            
            fig = figure('Name', ['Optimization Mux ' num2str(M)], 'Color', 'w', 'Position', [100, 100, 800, 600], 'Visible', fig_vis);
            
            subplot(2,1,1);
            plot(1:k, CNr(1:k), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
            hold on; 
            yline(best_cond, '--r', ['Min: ' num2str(best_cond, '%.2f')], 'LabelVerticalAlignment', 'bottom');
            set(gca, 'YScale', 'log', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k');
            box on; title(['Condition Number Evolution (M=' num2str(M) ')'], 'Color', 'k');
            ylabel('Cond(C)'); xlabel('Iterations');
            
            subplot(2,1,2);
            plot(1:k, sumden(1:k), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
            hold on;
            yline(best_density, '--r', ['Max: ' num2str(best_density, '%.4f')], 'LabelVerticalAlignment', 'bottom');
            set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k');
            box on; title('Density Evolution', 'Color', 'k');
            ylabel('Density'); xlabel('Iterations');
            
            if plot_mode == 1
                if ~exist(export_path, 'dir'), mkdir(export_path); end
                filename = fullfile(export_path, ['Opt_Mux_', num2str(M), '.png']);
                saveas(fig, filename);
                close(fig);
                fprintf('>> Plot M=%d exported to: %s\n', M, filename);
            end
        end

        C_optimized = best_C; % Save optimized matrix

        % Calculate optimized uniformity RMSE
        sum_x_opt = reshape(sum(C_optimized, [2, 3]), [N, 1]);
        sum_y_opt = reshape(sum(C_optimized, [1, 3]), [N, 1]);
        sum_z_opt = reshape(sum(C_optimized, [1, 2]), [N, 1]);
        errors_opt = [(sum_x_opt - M); (sum_y_opt - M); (sum_z_opt - M)];
        optimized_rmse = sqrt(mean(errors_opt.^2));


        %% --- UNIFORMIZATION PROCESS (POST-OPTIMIZATION) ---
        fprintf('>> Applying uniformization process to the best matrix found...\n');
        
        % Capture "Before" state for visualization
        C_before = best_C;

        % Apply uniformizer
        best_C = uniformizeMatrix(best_C, M);
        C_uniformized = best_C; % Save uniformized matrix
        final_cond = GetCond(best_C);

        % Calculate density after uniformization
        if R == 1
            [density_after, ~] = ComputeDensityIrregularSP(best_C);
        else
            [density_after, ~] = ComputeDensityRegularSP(best_C);
        end
        fprintf('>> Condition after uniformization: %.4f | Density: %.4f\n', final_cond, density_after);

        %% --- UNIFORMITY ANALYSIS ---
        % Sums by vectors (lines) for text report
        sum_lines_x = reshape(sum(best_C, [2, 3]), [N, 1]);
        sum_lines_y = reshape(sum(best_C, [1, 3]), [N, 1]);
        sum_lines_z = reshape(sum(best_C, [1, 2]), [N, 1]);
        
        non_uniform_x = sum(sum_lines_x ~= M);
        non_uniform_y = sum(sum_lines_y ~= M);
        non_uniform_z = sum(sum_lines_z ~= M);
        
        errors_final = [(sum_lines_x - M); (sum_lines_y - M); (sum_lines_z - M)];
        final_rmse = sqrt(mean(errors_final.^2));

        fprintf('>> Non-uniform vectors (R,C,L): %d, %d, %d | Uniformity RMSE: %.4f\n', ...
            non_uniform_x, non_uniform_y, non_uniform_z, final_rmse);

        if plot_mode > 0
            % --- Calculate global color limits for a fair comparison ---
            s_b_yz = squeeze(sum(C_before, 1));
            s_b_xz = squeeze(sum(C_before, 2));
            s_b_xy = squeeze(sum(C_before, 3));
            s_a_yz = squeeze(sum(best_C, 1));
            s_a_xz = squeeze(sum(best_C, 2));
            s_a_xy = squeeze(sum(best_C, 3));
            
            all_sums = [s_b_yz(:); s_b_xz(:); s_b_xy(:); s_a_yz(:); s_a_xz(:); s_a_xy(:)];
            color_lims = [min(all_sums), max(all_sums)];
            if color_lims(1) == color_lims(2), color_lims(2) = color_lims(1) + 1; end % Avoid error if all values are the same

            fig_vis = 'on';
            if plot_mode == 1, fig_vis = 'off'; end

            fig_uniformity = figure('Name', ['Projections Before and After Mux ' num2str(M)], 'Color', 'w', 'Position', [100, 100, 1200, 700], 'Visible', fig_vis);
            
            % --- BEFORE ---
            sgtitle('Uniformity Analysis', 'Color', 'k', 'FontSize', 16);
            
            ax1 = subplot(2,3,1);
            imagesc(s_b_yz); axis image; cb1 = colorbar;
            title('Before: YZ Projection', 'Color', 'k'); xlabel('Z'); ylabel('Y');
            set(ax1, 'XColor', 'k', 'YColor', 'k');
            set(cb1, 'Color', 'k');
            clim(ax1, color_lims);

            ax2 = subplot(2,3,2);
            imagesc(s_b_xz); axis image; cb2 = colorbar;
            title('Before: XZ Projection', 'Color', 'k'); xlabel('Z'); ylabel('X');
            set(ax2, 'XColor', 'k', 'YColor', 'k');
            set(cb2, 'Color', 'k');
            clim(ax2, color_lims);

            ax3 = subplot(2,3,3);
            imagesc(s_b_xy); axis image; cb3 = colorbar;
            title('Before: XY Projection', 'Color', 'k'); xlabel('Y'); ylabel('X');
            set(ax3, 'XColor', 'k', 'YColor', 'k');
            set(cb3, 'Color', 'k');
            clim(ax3, color_lims);

            % --- AFTER ---
            ax4 = subplot(2,3,4);
            imagesc(s_a_yz); axis image; cb4 = colorbar;
            title('After: YZ Projection', 'Color', 'k'); xlabel('Z'); ylabel('Y');
            set(ax4, 'XColor', 'k', 'YColor', 'k');
            set(cb4, 'Color', 'k');
            clim(ax4, color_lims);

            ax5 = subplot(2,3,5);
            imagesc(s_a_xz); axis image; cb5 = colorbar;
            title('After: XZ Projection', 'Color', 'k'); xlabel('Z'); ylabel('X');
            set(ax5, 'XColor', 'k', 'YColor', 'k');
            set(cb5, 'Color', 'k');
            clim(ax5, color_lims);

            ax6 = subplot(2,3,6);
            imagesc(s_a_xy); axis image; cb6 = colorbar;
            title('After: XY Projection', 'Color', 'k'); xlabel('Y'); ylabel('X');
            set(ax6, 'XColor', 'k', 'YColor', 'k');
            set(cb6, 'Color', 'k');
            clim(ax6, color_lims);

            if plot_mode == 1
                filename_uniformity = fullfile(export_path, ['Uniformity_Compare_Mux_', num2str(M), '.png']);
                saveas(fig_uniformity, filename_uniformity);
                close(fig_uniformity);
                fprintf('>> Uniformity comparison plot M=%d exported to: %s\n', M, filename_uniformity);
            end

            % --- Sphere Visualization Plot ---
            if show_sphere_plots
                % --- Custom Sphere Comparison Plot ---
                fig_comp = figure('Name', ['Sphere Compare Mux ' num2str(M)], 'Color', 'w', 'Position', [100, 100, 1800, 600], 'Visible', fig_vis);
                
                % Subplot 1: Initial State (Regular Spheres)
                subplot(1,3,1);
                showSphereCAReg(C_initial, 12);
                title(['Inicial (Esferas Regulares, M=' num2str(M) ')'], 'Color', 'k');
                
                % Subplot 2: Optimized State (Irregular Spheres)
                subplot(1,3,2);
                showSphereCA(C_optimized, 12, 1); % op=1 for irregular
                title(['Optimizado (Cond=' num2str(best_cond, '%.2f') ')'], 'Color', 'k');

                % Subplot 3: Uniformized State (Irregular Spheres)
                subplot(1,3,3);
                showSphereCA(C_uniformized, 12, 1); % op=1 for irregular
                title(['Uniformizado (Cond=' num2str(final_cond, '%.2f') ')'], 'Color', 'k');

                if plot_mode == 1
                    filename = fullfile(export_path, ['Sphere_Compare_Mux_', num2str(M), '.png']);
                    saveas(fig_comp, filename);
                    close(fig_comp);
                    fprintf('>> Gráfico de comparación de esferas M=%d exportado a: %s\n', M, filename);
                end
            end
        end
    else
        C_optimized = best_C;   % For M=1, optimized is same as initial
        C_uniformized = best_C; % For M=1, uniformized is same as initial
        optimized_rmse = initial_rmse; % For M=1
        final_rmse = initial_rmse;     % For M=1
        fprintf('M=1 detected: Skipping optimization and uniformization.\n');
        final_cond = best_cond;

        if plot_mode > 0 && show_sphere_plots
            fig_vis = 'on';
            if plot_mode == 1, fig_vis = 'off'; end
            
            % --- Custom Sphere Comparison Plot ---
            fig_comp = figure('Name', ['Sphere Compare Mux ' num2str(M)], 'Color', 'w', 'Position', [100, 100, 1800, 600], 'Visible', fig_vis);
            
            % Subplot 1: Initial State (Regular Spheres)
            subplot(1,3,1);
            showSphereCAReg(C_initial, 12);
            title(['Inicial (Esferas Regulares, M=' num2str(M) ')'], 'Color', 'k');
            
            % Subplot 2: Optimized State (Irregular Spheres)
            subplot(1,3,2);
            showSphereCA(C_optimized, 12, 1); % op=1 for irregular
            title(['Optimizado (Cond=' num2str(best_cond, '%.2f') ')'], 'Color', 'k');

            % Subplot 3: Uniformized State (Irregular Spheres)
            subplot(1,3,3);
            showSphereCA(C_uniformized, 12, 1); % op=1 for irregular
            title(['Uniformizado (Cond=' num2str(final_cond, '%.2f') ')'], 'Color', 'k');

            if plot_mode == 1
                filename = fullfile(export_path, ['Sphere_Compare_Mux_', num2str(M), '.png']);
                saveas(fig_comp, filename);
                close(fig_comp);
                fprintf('>> Gráfico de comparación de esferas M=%d exportado a: %s\n', M, filename);
            end
        end
    end
    %% --- Final Matrix Processing ---
    N1 = ceil(N/NF);
    II = ones(N1,N1);
    T = zeros(N, N, NF);
    for i_f = 1:NF
        temp_kron = kron(II, best_C(:,:,i_f));
        T(:,:,i_f) = temp_kron(1:N, 1:N);
    end
    C_final = T;

    fprintf('-------------------- Resumen M=%d --------------------\n', M);
    fprintf('  Condición: Inicial: %.2f | Optimizada: %.2f | Uniforme: %.4f\n', initial_cond, best_cond, final_cond);
    fprintf('  Densidad:  Inicial: %.4f | Optimizada: %.4f | Uniforme: %.4f\n', initial_density, best_density, density_after);
    fprintf('  RMSE Unif: Inicial: %.4f | Optimizada: %.4f | Uniforme: %.4f\n\n', initial_rmse, optimized_rmse, final_rmse);
    
    % Escribir el mismo resumen en el archivo de log
    fprintf(fileID, '-------------------- Resumen M=%d --------------------\n', M);
    fprintf(fileID, '  Condición: Inicial: %.2f | Optimizada: %.2f | Uniforme: %.4f\n', initial_cond, best_cond, final_cond);
    fprintf(fileID, '  Densidad:  Inicial: %.4f | Optimizada: %.4f | Uniforme: %.4f\n', initial_density, best_density, density_after);
    fprintf(fileID, '  RMSE Unif: Inicial: %.4f | Optimizada: %.4f | Uniforme: %.4f\n\n', initial_rmse, optimized_rmse, final_rmse);
    
    % --- Guardado de Resultados ---
    if sf==1
        % Crear structs para un guardado más claro
        initial_pattern.matrix = C_initial;
        initial_pattern.condition_number = initial_cond;
        initial_pattern.density = initial_density;
        initial_pattern.uniformity_rmse = initial_rmse;
        
        optimized_pattern.matrix = C_optimized;
        optimized_pattern.condition_number = best_cond;
        optimized_pattern.density = best_density;
        optimized_pattern.uniformity_rmse = optimized_rmse;
        
        uniformized_pattern.matrix = C_uniformized;
        uniformized_pattern.condition_number = final_cond;
        uniformized_pattern.density = density_after;
        uniformized_pattern.uniformity_rmse = final_rmse;

        folderName = "./Patrones/SSPv4";
        if ~exist(folderName, 'dir'), mkdir(folderName); end
        save(fullfile(folderName, ['Bands_', num2str(NF), '_Mux_', num2str(M), '.mat']), ...
            'initial_pattern', 'optimized_pattern', 'uniformized_pattern');
    end

    %% --- Store metrics for summary plot ---
    summary_initial_cond(m_idx) = initial_cond;
    summary_optimized_cond(m_idx) = best_cond;
    summary_uniform_cond(m_idx) = final_cond;
    summary_initial_density(m_idx) = initial_density;
    summary_optimized_density(m_idx) = best_density;
    summary_uniform_density(m_idx) = density_after;
    summary_initial_rmse(m_idx) = initial_rmse;
    summary_optimized_rmse(m_idx) = optimized_rmse;
    summary_uniform_rmse(m_idx) = final_rmse;
end

% --- Close log file ---
fclose(fileID);
fprintf('>> Execution summary saved to: %s\n', log_filename);

%% --- Final Summary Plot (Comparison between M) ---
if plot_mode > 0 && length(mux_values_to_run) > 1
    fig_vis = 'on';
    if plot_mode == 1, fig_vis = 'off'; end

    summary_fig = figure('Name', 'Metrics Summary vs. M', 'Color', 'w', 'Position', [100, 100, 900, 900], 'Visible', fig_vis);
    
    % Subplot for Condition Number
    subplot(3,1,1);
    
    % --- Visualization logic for very large and Inf values ---
    plot_initial = summary_initial_cond;
    plot_optimized = summary_optimized_cond;
    plot_uniform = summary_uniform_cond;

    % Identify values that should go to the top of the plot
    is_initial_huge = isinf(plot_initial) | (plot_initial > 10000 * max(plot_optimized, plot_uniform));
    is_optimized_inf = isinf(plot_optimized);
    is_uniform_inf = isinf(plot_uniform);

    % Determine Y-axis limits based only on finite and "reasonable" values
    sane_vals = [plot_initial(~is_initial_huge), plot_optimized(~is_optimized_inf), plot_uniform(~is_uniform_inf)];
    if isempty(sane_vals) || all(isnan(sane_vals)), sane_vals = 1; end
    min_sane_y = min(sane_vals(sane_vals>0));
    if isempty(min_sane_y), min_sane_y = 1; end
    max_sane_y = max(sane_vals);
    y_limits = [max(1e-1, min_sane_y / 10), min(1e12, max_sane_y * 10)];
    if y_limits(1) >= y_limits(2), y_limits(1) = y_limits(2)/100; end

    % Replace huge/inf values with a value above the limit so the line is drawn upwards
    replacement_val = y_limits(2) * 100; % Increased to ensure it's well outside
    plot_initial(is_initial_huge) = replacement_val;
    plot_optimized(is_optimized_inf) = replacement_val;
    plot_uniform(is_uniform_inf) = replacement_val;

    % Draw the lines
    p1 = plot(mux_values_to_run, plot_initial, 'o-', 'LineWidth', 2, 'DisplayName', 'Initial');
    hold on;
    p2 = plot(mux_values_to_run, plot_optimized, 's-', 'LineWidth', 2, 'DisplayName', 'Optimized');
    p3 = plot(mux_values_to_run, plot_uniform, '^-', 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880], 'DisplayName', 'Uniform');
    hold off;

    % Set scale and limits to "clip" the view
    set(gca, 'YScale', 'log', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    ylim(y_limits);

    title('Condition Number Evolution vs. M', 'Color', 'k');
    xlabel('M Value (Mux)');
    ylabel('Condition Number');
    legend('show', 'Location', 'northwest', 'Color', 'white', 'TextColor', 'black');
    grid on;

    % Subplot for Density
    subplot(3,1,2);
    plot(mux_values_to_run, summary_initial_density, 'o-', 'LineWidth', 2, 'DisplayName', 'Initial');
    hold on;
    plot(mux_values_to_run, summary_optimized_density, 's-', 'LineWidth', 2, 'DisplayName', 'Optimized');
    plot(mux_values_to_run, summary_uniform_density, '^-', 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880], 'DisplayName', 'Uniform');
    hold off;
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    title('Density Evolution vs. M', 'Color', 'k');
    xlabel('M Value (Mux)');
    ylabel('Density');
    legend('show', 'Location', 'best', 'Color', 'white', 'TextColor', 'black');
    grid on;

    % Subplot for Uniformity RMSE
    subplot(3,1,3);
    plot(mux_values_to_run, summary_initial_rmse, 'o-', 'LineWidth', 2, 'DisplayName', 'Initial');
    hold on;
    plot(mux_values_to_run, summary_optimized_rmse, 's-', 'LineWidth', 2, 'DisplayName', 'Optimized');
    plot(mux_values_to_run, summary_uniform_rmse, '^-', 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880], 'DisplayName', 'Uniform');
    hold off;
    set(gca, 'YScale', 'log', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    title('Uniformity RMSE Evolution vs. M', 'Color', 'k');
    xlabel('M Value (Mux)');
    ylabel('Uniformity RMSE');
    legend('show', 'Location', 'best', 'Color', 'white', 'TextColor', 'black');
    grid on;

    if plot_mode == 1
        filename_summary = fullfile(export_path, 'Summary_Metrics_vs_M.png');
        saveas(summary_fig, filename_summary);
        close(summary_fig);
        fprintf('>> Summary plot exported to: %s\n', filename_summary);
    end
end

function C_out = uniformizeMatrix(C_in, M)
    % uniformizeMatrix: Adjusts a 3D matrix so that the sum of each line is M.
    % The strategy is to accept/reject random swaps based on RMSE improvement.
    C_out = C_in;
    [N, ~, ~] = size(C_out);
    max_iter = 200 * N^2; % Increased number of attempts to balance

    for iter = 1:max_iter
        % 1. Calculate line sums and difference from target M
        sum_x = reshape(sum(C_out, [2, 3]), [N, 1]); % Sum by rows (X-axis)
        sum_y = reshape(sum(C_out, [1, 3]), [N, 1]); % Sum by columns (Y-axis)
        sum_z = reshape(sum(C_out, [1, 2]), [N, 1]); % Sum by layers (Z-axis)

        diff_x = sum_x - M;
        diff_y = sum_y - M;
        diff_z = sum_z - M;

        % Calculate current RMSE for comparison
        current_rmse = sqrt(mean([diff_x.^2; diff_y.^2; diff_z.^2]));

        % 2. If everything is balanced, the job is done.
        if ~any(diff_x) && ~any(diff_y) && ~any(diff_z)
            return;
        end

        % 3. Create a mask to block "perfect" points.
        % A point (i,j,k) is blocked if its row i, AND its column j, AND its layer k already sum to M.
        is_row_ok = (diff_x == 0);
        is_col_ok = (diff_y == 0);
        is_layer_ok = (diff_z == 0);

        % Use implicit expansion (broadcasting) to create the 3D mask
        LockMask = reshape(is_row_ok, N, 1, 1) & reshape(is_col_ok, 1, N, 1) & reshape(is_layer_ok, 1, 1, N);

        % 4. Identify '1's and '0's that can be moved.
        % These are the points that are NOT in the lock mask.
        movable_ones_indices = find(C_out & ~LockMask);
        fillable_zeros_indices = find(~C_out & ~LockMask);

        % 5. If no moves are possible, the algorithm is stuck.
        if isempty(movable_ones_indices) || isempty(fillable_zeros_indices)
            warning('Uniformizer: No valid moves found. Matrix may not be perfectly uniform.');
            return;
        end

        % 6. Perform a random swap, protecting against zero sums.
        [ones_i, ones_j, ones_k] = ind2sub(size(C_out), movable_ones_indices);
        % A '1' is critical if it is the only one in its row, column, or layer.
        is_critical = (sum_x(ones_i) == 1) | (sum_y(ones_j) == 1) | (sum_z(ones_k) == 1);
        safe_movable_ones = movable_ones_indices(~is_critical);
        
        indices_to_use_for_one = safe_movable_ones;
        if isempty(indices_to_use_for_one)
            indices_to_use_for_one = movable_ones_indices; % Fallback if all are critical
        end

        rand_one_idx = randi(length(indices_to_use_for_one));
        rand_zero_idx = randi(length(fillable_zeros_indices));

        one_to_move_flat_idx = indices_to_use_for_one(rand_one_idx);
        zero_to_fill_flat_idx = fillable_zeros_indices(rand_zero_idx);

        % 7. Perform the swap on a temporary copy and evaluate if it improves the RMSE.
        C_temp = C_out;
        C_temp(one_to_move_flat_idx) = 0;
        C_temp(zero_to_fill_flat_idx) = 1;

        % Recalculate RMSE for the temporary matrix
        sum_x_new = reshape(sum(C_temp, [2, 3]), [N, 1]);
        sum_y_new = reshape(sum(C_temp, [1, 3]), [N, 1]);
        sum_z_new = reshape(sum(C_temp, [1, 2]), [N, 1]);
        errors_new = [(sum_x_new - M); (sum_y_new - M); (sum_z_new - M)];
        new_rmse = sqrt(mean(errors_new.^2));

        % 8. Accept the change only if it improves (or maintains) the RMSE.
        if new_rmse <= current_rmse
            C_out = C_temp;
        end
    end
    warning('Uniformizer: Maximum iterations reached. Matrix may not be perfectly uniform.');
end