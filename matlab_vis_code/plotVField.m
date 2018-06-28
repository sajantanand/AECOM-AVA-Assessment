function [] = plotVField(xPos,yPos,zPos,xVel,yVel,zVel,velMag)
%% params
lw = 1.5; % linewidth of arrows
hs = 2; % max head size for arrows

%% plotting
hold on
len =  size(xPos,1);
Q=quiver3(xPos,yPos,zPos,  xVel,yVel,zVel, 'LineWidth',lw,'maxheadsize',hs);
%S=scatter3(xPos,yPos,repmat(2,len,1),50,velMag,'filled');

%% color based on total magnitude (note to Matlab: please make such things simpler to do)
%// Get the current colormap
currentColormap = colormap(gca,jet);
ax = gca;
set(ax,'Color',[.4 .4 .4]);

%// Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(velMag, size(currentColormap, 1));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(Q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(Q.Tail,'ColorBinding', 'interpolated','ColorData', reshape(cmap(1:2,:,:), [], 4).');
end