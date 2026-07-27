% -------------------------------------------------------------------------
% Script: generate_comparison_plots.m
% Autor: Gemini Code Assist
% Date:   July 2026
%
% Description:
% This script generates a comparative plot of three sphere packing states for
% specific M values (8 and 15), loading data from specified folders:
% 1. Random State: Loaded from the 'nUnif' directory.
% 2. Initial State: Loaded from the 'initial_pattern' struct in the SSPv4 file.
% 3. Optimized State: Loaded from the 'NewDesignSSP' directory.
%
% The script replicates the 'Sphere Compare Mux' plots from GenCSSPv4_NF31.
% -------------------------------------------------------------------------

% Añadir las rutas necesarias para acceder a las funciones de utilidad y ploteo
addpath(genpath('./src'));
addpath(genpath('./utils'));
addpath(genpath('./Method'));
addpath(genpath('./SubKernels'));
addpath(genpath('./PlotSpheres'));
clc; clear; close all;

%% --- Configuration ---
NF = 31; % Bands (Z-depth), must match the saved patterns
mux_values_to_run = [8, 15];
export_path = './comparison_graphs'; % Directory to save the plots

% --- Visualization Settings ---
fig_vis = 'on'; % 'on' to show plots, 'off' for headless execution
save_plots = true;

%% --- Main Loop ---
for m_idx = 1:length(mux_values_to_run)
    M = mux_values_to_run(m_idx);
    fprintf('--- Processing M = %d ---\n', M);

    % --- File Paths ---
    random_file_path = fullfile('C:\Users\ale_m\Escritorio\ShuffledSpherePacking\Patrones\nUnif', ...
                                   ['Bands_', num2str(NF), '_Mux_', num2str(M), '.mat']);

    sspv4_file_path = fullfile('C:\Users\ale_m\Escritorio\ShuffledSpherePacking\Patrones\SSPv4', ...
                                   ['Bands_', num2str(NF), '_Mux_', num2str(M), '.mat']);

    optimized_file_path = fullfile('C:\Users\ale_m\Escritorio\ShuffledSpherePacking\Patrones\NewDesignSSP', ...
                                   ['Bands_', num2str(NF), '_Mux_', num2str(M), '.mat']);

    % --- 1. Load Random State (from nUnif) ---
    fprintf('Loading Random state from: %s\n', random_file_path);
    if exist(random_file_path, 'file')
        loaded_random = load(random_file_path);
        % Robustly load the matrix from the file
        if isfield(loaded_random, 'C')
             C_random = loaded_random.C;
        elseif isfield(loaded_random, 'C_final')
             C_random = loaded_random.C_final;
        else
            fprintf('ERROR: Could not find a recognized matrix variable in the Random state file.\n');
            continue;
        end
        fprintf('  -> Random state matrix loaded from nUnif.\n');
    else
        fprintf('ERROR: Random state file not found. Skipping M=%d.\n', M);
        continue;
    end

    % --- 2. Load Initial and Optimized States ---
    fprintf('Loading Initial state from: %s\n', sspv4_file_path);
    if exist(sspv4_file_path, 'file')
        loaded_sspv4 = load(sspv4_file_path);
        
        % --- Process Initial State (from SSPv4) ---
        if isfield(loaded_sspv4, 'initial_pattern')
            C_initial_sspv4 = loaded_sspv4.initial_pattern.matrix;
            initial_cond = loaded_sspv4.initial_pattern.condition_number;
            [initial_density, initial_diameter] = ComputeDensityIrregularSP(C_initial_sspv4);
            fprintf('  -> Initial (SSPv4) pattern loaded.\n');
        else
            fprintf('ERROR: "initial_pattern" struct not found in SSPv4 file.\n');
            continue; % Skip to the next M
        end
    else
        fprintf('ERROR: SSPv4 state file not found for Initial pattern. Skipping M=%d.\n', M);
        continue;
    end

    fprintf('Loading Optimized state from: %s\n', optimized_file_path);
    if exist(optimized_file_path, 'file')
        loaded_optimized = load(optimized_file_path);
        % Robustly load the matrix from the file
        if isfield(loaded_optimized, 'C')
            C_optimized = loaded_optimized.C;
        elseif isfield(loaded_optimized, 'C_final')
            C_optimized = loaded_optimized.C_final;
        elseif isfield(loaded_optimized, 'optimized_pattern')
            C_optimized = loaded_optimized.optimized_pattern.matrix;
        else
            fprintf('ERROR: Could not find a recognized matrix variable in the Optimized state file.\n');
            continue; % Skip to the next M
        end
        fprintf('  -> Optimized state matrix loaded from NewDesignSSP.\n');
    else
        fprintf('ERROR: Optimized state file not found. Skipping M=%d.\n', M);
        continue;
    end

    % --- 3. Calculate Metrics for Titles ---
    fprintf('Calculating metrics for plots...\n');
    % Metrics for Random state
    random_cond = GetCond(C_random);
    [random_density, random_diameter] = ComputeDensityIrregularSP(C_random);

    % Metrics for Optimized state
    optimized_cond = GetCond(C_optimized);
    [optimized_density, optimized_diameter] = ComputeDensityIrregularSP(C_optimized);

    % --- 4. Generate Comparative Plot ---
    fprintf('Generating sphere comparison plot...\n');
    fig_comp = figure('Name', ['Sphere Comparison Mux ' num2str(M)], 'Position', [100, 100, 1800, 600], 'Visible', fig_vis);

    % Subplot 1: Random State (from nUnif)
    subplot(1,3,1);
    showSphereCA(C_random, 12,2);
    title({'Random State (from nUnif)', ...
           sprintf('Mux=%d, Cond=%.2f, Dens=%.4f, Diam=%.2f', M, random_cond, random_density, random_diameter)}, ...
          'Color', 'k');

    % Subplot 2: Initial State (Loaded from SSPv4)
    subplot(1,3,2);
    showSphereCA(C_initial_sspv4, 12, 1); % op=1 for irregular
    title({'Initial State (from SSPv4)', ...
           sprintf('Mux=%d, Cond=%.2f, Dens=%.4f, Diam=%.2f', M, initial_cond, initial_density, initial_diameter)}, ...
          'Color', 'k');

    % Subplot 3: Optimized State (Loaded from SSPv4)
    subplot(1,3,3);
    showSphereCA(C_optimized, 12, 1); % op=1 for irregular
    title({'Optimized State (from NewDesignSSP)', ...
           sprintf('Mux=%d, Cond=%.2f, Dens=%.4f, Diam=%.2f', M, optimized_cond, optimized_density, optimized_diameter)}, ...
          'Color', 'k');

    % --- 5. Apply Final Styling (White Background, Black Text) ---
    % set(fig_comp, 'Color', 'w'); % Set figure background to white
    % allAxes = findall(fig_comp, 'type', 'axes');
    % for ax = allAxes'
    %     ax.Color = 'w';    % Set plot background to white
    %     ax.XColor = 'k';   % Set X-axis ticks and label color to black
    %     ax.YColor = 'k';   % Set Y-axis ticks and label color to black
    %     ax.ZColor = 'k';   % Set Z-axis ticks and label color to black
    %     ax.GridColor = 'k';% Set grid color to black
    %     ax.GridAlpha = 0.25;
    % end

    if save_plots
        if ~exist(export_path, 'dir'), mkdir(export_path); end
        filename = fullfile(export_path, ['Sphere_Comparison_Mux_', num2str(M), '.png']);
        saveas(fig_comp, filename);
        fprintf('>> Plot saved to: %s\n', filename);
        if strcmpi(fig_vis, 'off')
            close(fig_comp);
        end
    end
    fprintf('--- M = %d Finished ---\n\n', M);
end

fprintf('All processes have finished.\n');