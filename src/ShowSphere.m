function [] = ShowSphere(G,radius,fontsize,NF)
[r,c,z1] = find(G);
C =colormap(jet(NF));
%markerColors = jet(numel(G));
for i=1:nnz(G)
    % Make unit sphere
    [x,y,z] = sphere(20);
    
    % Scale to desire radius.
    x = x * radius;
    y = y * radius;
    z = z * radius;
    % Translate sphere to new location.
    %offset = r(i), c(i), z1(i);
    % Plot as surface.
    
    surf(x+r(i),y+c(i),z+z1(i),'FaceColor', C(z1(:),:),'LineStyle',":");%, ...
       % 'LineWidth',0.1,'LineStyle',":",'EdgeColor',[0.3 0.3 0.3]); 
    pbaspect([1 1 1])
    % Label axes.
    hold on
end
end