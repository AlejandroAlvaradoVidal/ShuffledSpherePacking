% -------------------------------------------------------------------------
% Script: plot_loaded_matrices.m
% Autor: Gemini Code Assist
% Date:   July 2026
%
% Description:
%   Genera un gráfico comparativo para tres matrices 3D pre-cargadas:
%   1. C_rand: Matriz del estado aleatorio.
%   2. C_init: Matriz del estado inicial (del SSPv4).
%   3. C_optimized: Matriz del estado optimizado (del SSPv4).
%
% INSTRUCCIONES DE USO:
%   1. Carga tus tres matrices (C_rand, C_init, C_optimized) en el
%      espacio de trabajo de MATLAB.
%   2. Ajusta el valor de 'M' en la sección de configuración de abajo.
%   3. Ejecuta este script desde la ventana de comandos de MATLAB.
% -------------------------------------------------------------------------

%% --- Configuración del Usuario ---
% IMPORTANTE: Ajusta el valor de M para que coincida con tus datos cargados.
M = 8; % <--- CAMBIA ESTE VALOR SI ES NECESARIO

% --- Configuración de Visualización y Guardado ---
export_path = './comparison_graphs_from_workspace'; % Carpeta para guardar el gráfico
fig_vis = 'on'; % 'on' para mostrar el gráfico, 'off' para no mostrarlo
save_plots = true;

%% --- Verificación y Preparación ---
% Añadir las rutas a las funciones necesarias
addpath(genpath('./src'));
addpath(genpath('./utils'));
addpath(genpath('./Method'));
addpath(genpath('./SubKernels'));
addpath(genpath('./PlotSpheres'));

% Comprobar si las variables existen en el espacio de trabajo
if ~exist('C_rand', 'var') || ~exist('C_init', 'var') || ~exist('C_optimized', 'var')
    error('Error: Por favor, carga las matrices C_rand, C_init y C_optimized en el espacio de trabajo antes de ejecutar el script.');
end

fprintf('--- Procesando M = %d desde las variables del espacio de trabajo ---\n', M);

%% --- 1. Cálculo de Métricas para los Títulos ---
fprintf('Calculando métricas para los gráficos...\n');

% Métricas para el estado Aleatorio
random_cond = GetCond(C_rand);
[random_density, random_diameter] = ComputeDensityIrregularSP(C_rand);

% Métricas para el estado Inicial
initial_cond = GetCond(C_init);
[initial_density, initial_diameter] = ComputeDensityIrregularSP(C_init);

% Métricas para el estado Optimizado
optimized_cond = GetCond(C_optimized);
[optimized_density, optimized_diameter] = ComputeDensityIrregularSP(C_optimized);

%% --- 2. Generación del Gráfico Comparativo ---
fprintf('Generando el gráfico de comparación de esferas...\n');
fig_comp = figure('Name', ['Sphere Comparison Mux ' num2str(M)], 'Position', [100, 100, 1800, 600], 'Visible', fig_vis);

% Subplot 1: Estado Aleatorio
subplot(1,3,1);
showSphereCAReg(C_rand, 12);
title_str_1 = sprintf('Random State (Mux=%d, Cond=%.2f, Dens=%.4f, Diam=%.2f)', M, random_cond, random_density, random_diameter);
title(title_str_1, 'Color', 'k');

% Subplot 2: Estado Inicial
subplot(1,3,2);
showSphereCA(C_init, 12, 1); % op=1 para irregular
title_str_2 = sprintf('Initial State (Mux=%d, Cond=%.2f, Dens=%.4f, Diam=%.2f)', M, initial_cond, initial_density, initial_diameter);
title(title_str_2, 'Color', 'k');

% Subplot 3: Estado Optimizado
subplot(1,3,3);
showSphereCA(C_optimized, 12, 1); % op=1 para irregular
title_str_3 = sprintf('Optimized State (Mux=%d, Cond=%.2f, Dens=%.4f, Diam=%.2f)', M, optimized_cond, optimized_density, optimized_diameter);
title(title_str_3, 'Color', 'k');

%% --- 3. Aplicar Estilo Final (Fondo Blanco, Texto Negro) ---
set(fig_comp, 'Color', 'w');
allAxes = findall(fig_comp, 'type', 'axes');
for ax = allAxes'
    ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k'; ax.ZColor = 'k';
    ax.GridColor = 'k'; ax.GridAlpha = 0.25;
end

%% --- 4. Guardar el Gráfico ---
if save_plots
    if ~exist(export_path, 'dir'), mkdir(export_path); end
    filename = fullfile(export_path, ['Workspace_Comparison_Mux_', num2str(M), '.png']);
    saveas(fig_comp, filename);
    fprintf('>> Gráfico guardado en: %s\n', filename);
    if strcmpi(fig_vis, 'off'), close(fig_comp); end
end

fprintf('--- Generación de gráfico finalizada ---\n\n');