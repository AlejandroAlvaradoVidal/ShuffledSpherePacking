function [] = showSphereGeneral(G,fontsize)
[r,c,z1] = find(G);
X = [r c z1]'; % generator matrix
D = pdist(X');
D1 = squareform(D);
%D1= reshape(D,size(G,1)^2,[]);
R = zeros(numel(G),1);
for i=1:numel(G)
     temp = D1(i,:);
     [x2,y2,z2]=find(temp);
     R(i)= min(z2)/2;
end
%markerColors = jet(numel(G));
for i=1:numel(G)
    % Make unit sphere
    [x,y,z] = sphere(20);
    
    % Scale to desire radius.
    x = x * R(i);
    y = y * R(i);
    z = z * R(i);
    % Translate sphere to new location.
    %offset = r(i), c(i), z1(i);
    % Plot as surface.
    C =colormap(jet(max(z1)));
    surf(x+r(i),y+c(i),z+z1(i),'FaceColor', C(z1(i),:)) 
    pbaspect([1 1 1])
    % Label axes.
    hold on
end
axis equal
end