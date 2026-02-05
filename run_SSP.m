function [C_final, best_cond] = optimizeCodedAperture(NF, N, M, options)
    % optimizeCodedAperture: Optimiza la apertura codificada para un M específico.
    
    arguments
        NF (1,1) double
        N (1,1) double
        M (1,1) double
        options.R (1,1) double = 1
        options.max_iterations (1,1) double = 2000
        options.rel_tol (1,1) double = 1e-5
        options.patience (1,1) double = 60
        options.burn_in (1,1) double = 300
        options.plot_mode (1,1) double = 1
        options.export_path char = './graphs'
        options.save_file (1,1) logical = false
    end

    % --- Inicialización ---
    [G] = RunEQ10MUX4(NF, M);
    [C, ~] = generateCodedAperture(G, NF, M, NF); 
    C(C > 1) = 1;
    
    Cond = GetCond(C);
    best_cond = Cond;
    best_C = C;
    
    if M <= 1
        fprintf('M=1 detectado: Saltando optimización.\n');
    else
        % Históricos
        CNr = nan(1, options.max_iterations); 
        sumden = nan(1, options.max_iterations);
        no_improvement_count = 0;
        
        fprintf('Iniciando Optimización M=%d | Cond Inicial: %.2f\n', M, Cond);
        
        for k = 1:options.max_iterations
            [list, diameters] = uniqueDistanceIrregularSP(C);
            
            % Extraer coordenadas de elementos activos
            X = [];
            for i_layer = 1:NF
                tp = C(:,:,i_layer) * i_layer;
                [x_idx, y_idx, z_idx] = find(tp);
                X = [X; [x_idx, y_idx, z_idx]];
            end
            
            positions = find(diameters == list(1));
            val = randi([1, length(positions)], 1, 1);
            [xyzt] = X(positions(val), :);
            i = xyzt(1); ii = xyzt(2); j = xyzt(3);
            
            C_next = C;
            if rand() < 0.5
                % INTERCAMBIO PLANO XZ
                otherLayer = randi([1, NF]);
                temp = C(i, :, j); 
                C(i, :, j) = C(i, :, otherLayer);
                C(i, :, otherLayer) = temp;
            else
                % INTERCAMBIO PLANO YZ
                otherLayer = randi([1, NF]);
                temp = C(:, ii, j);
                C(:, ii, j) = C(:, ii, otherLayer);
                C(:, ii, otherLayer) = temp;
            end
            C_next(C_next > 1) = 1;
            
            new_cond = GetCond(C_next);
            new_density = calculateDensity(C_next, options.R);

            if new_cond < best_cond
                improvement = (best_cond - new_cond) / best_cond;
                if improvement > options.rel_tol
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
                fprintf('M=%d | Iter: %d | Best: %.2f | Status: %d/%d\n', ...
                    M, k, best_cond, no_improvement_count, options.patience);
            end
            
            if k > options.burn_in && no_improvement_count >= options.patience
                break;
            end
        end
        
        % --- Visualización ---
        if options.plot_mode > 0
            renderPlots(CNr, sumden, k, M, best_cond, options);
        end
    end

    % --- Procesamiento Final ---
    N1 = ceil(N/NF);
    II = ones(N1,N1);
    C_final = zeros(N, N, NF);
    for i_f = 1:NF
        temp_kron = kron(II, best_C(:,:,i_f));
        C_final(:,:,i_f) = temp_kron(1:N, 1:N);
    end

    if options.save_file
        saveResults(C_final, NF, M);
    end
end

%% --- Funciones Auxiliares Internas ---

function d = calculateDensity(C, R)
    if R == 1
        [d, ~] = ComputeDensityIrregularSP(C);
    else
        [d, ~] = ComputeDensityRegularSP(C);
    end
end

function renderPlots(CNr, sumden, k, M, best_cond, options)
    fig_vis = 'on';
    if options.plot_mode == 1, fig_vis = 'off'; end
    
    fig = figure('Name', ['Opt Mux ' num2str(M)], 'Color', 'w', 'Visible', fig_vis);
    
    % Subplot 1: Condition Number
    subplot(2,1,1);
    plot(1:k, CNr(1:k), 'Color', [0 0.447 0.741], 'LineWidth', 1.5);
    hold on;
    yline(best_cond, '--r', sprintf('Min: %.2f', best_cond));
    set(gca, 'YScale', 'log', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    title(['Evolución Condition Number (M=' num2str(M) ')'], 'Color', 'k');
    grid on;

    % Subplot 2: Density
    subplot(2,1,2);
    plot(1:k, sumden(1:k), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    title('Evolución de la Densidad', 'Color', 'k');
    xlabel('Iteraciones'); grid on;

    if options.plot_mode == 1
        if ~exist(options.export_path, 'dir'), mkdir(options.export_path); end
        saveas(fig, fullfile(options.export_path, ['Opt_Mux_', num2str(M), '.png']));
        close(fig);
    end
end

function saveResults(C, NF, M)
    folderName = "./Patrones/SSPx2";
    if ~exist(folderName, 'dir'), mkdir(folderName); end
    save(fullfile(folderName, sprintf('Bands_%d_Mux_%d.mat', NF, M)), "C");
end