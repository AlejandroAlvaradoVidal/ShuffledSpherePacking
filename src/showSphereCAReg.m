function [] = showSphereCAReg(C,fontsize)

[list,xyz]= ComputeDistanceSP(C);
%R=mean(list)*ones(size(list));
R=mean(list).*ones(size(list));

cmap=colormap(jet(size(C,3)));
for i=1:nnz(C)
    % Make unit sphere
    [x,y,z] = sphere(10);
    
    % Scale to desire radius.
    x = x * R(i);
    y = y * R(i);
    z = z * R(i);
    % Translate sphere to new location.
    %offset = r(i), c(i), z1(i);
    % Plot as surface.
    %C =colormap(jet(max(z1)));
    surf(x+xyz(i,1),y+xyz(i,2),z+xyz(i,3),'FaceColor', cmap(xyz(i,3),:),'EdgeColor',[0.5 0.5 0.5],'LineStyle','-','LineWidth',0.5,'EdgeAlpha',0.2) 
    pbaspect([1 1 1])
    xlim([0,size(C,3)+1])
    ylim([0,size(C,3)+1])
    zlim([0,size(C,3)+1])
    % Label axes.
    hold on
end

% 1. Forzar el fondo de la ventana a blanco
fig = gcf;
fig.Color = 'w'; 

% 2. Encontrar TODOS los subplots/ejes de la figura actual
allAxes = findall(fig, 'type', 'axes');

% 3. Aplicar los colores claros a todos los ejes simultáneamente
set(allAxes, 'Color', 'w');       % Fondo del cubo 3D en blanco
set(allAxes, 'XColor', 'k');      % Texto y líneas del eje X en negro
set(allAxes, 'YColor', 'k');      % Texto y líneas del eje Y en negro
set(allAxes, 'ZColor', 'k');      % Texto y líneas del eje Z en negro

% (Opcional) Oscurecer la cuadrícula para que resalte bien sobre el blanco
set(allAxes, 'GridColor', 'k');
set(allAxes, 'GridAlpha', 0.15);  % Transparencia de la cuadrícula
axis equal
hold off
end