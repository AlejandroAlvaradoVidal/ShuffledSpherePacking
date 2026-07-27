function [] = showSphereComparison(C_initial, C_optimized, C_uniformized, M, fig_vis, export_path)
%showSphereComparison Creates and saves a figure with 3D sphere plots.
%   Generates a side-by-side comparison of the initial, optimized, and
%   uniformized matrices.

    fig = figure('Name', ['Sphere Visualization M=' num2str(M)], 'Color', 'w', 'Position', [50, 50, 1500, 500], 'Visible', fig_vis);

    
    
    % Subplot 1: Initial
    subplot(1,3,1);
    showSphereCA(C_initial, 10, 1);
    title(['Initial (M=' num2str(M) ')'], 'Color', 'k');

    % Subplot 2: Optimized
    subplot(1,3,2);
    showSphereCA(C_optimized, 10, 1);
    title(['Optimized (M=' num2str(M) ')'], 'Color', 'k');

    % Subplot 3: Uniformized
    subplot(1,3,3);
    showSphereCA(C_uniformized, 10, 1);
    title(['Uniformized (M=' num2str(M) ')'], 'Color', 'k');

    % Export if in headless mode
    if strcmp(fig_vis, 'off')
        if ~exist(export_path, 'dir'), mkdir(export_path); end
        filename = fullfile(export_path, ['Sphere_Compare_Mux_', num2str(M), '.png']);
        saveas(fig, filename);
        close(fig);
        fprintf('>> Sphere comparison plot M=%d exported to: %s\n', M, filename);
    end
end