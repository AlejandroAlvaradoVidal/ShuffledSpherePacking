% main_runner.m
addpath(genpath('./src'));
addpath(genpath('./utils'));
% ... agregar el resto de tus paths ...

% Parámetros base
NF = 31;
N = 31;
mux_values = [1 : 16]; % Puedes iterar sobre varios valores fácilmente

for M = mux_values
    [C_final, cond] = optimizeCodedAperture(NF, N, M, ...
        'max_iterations', 2000, ...
        'patience', 60, ...
        'plot_mode', 1, ...
        'save_file', true);
    
    fprintf('Finalizado M=%d con Cond: %.4f\n', M, cond);
end