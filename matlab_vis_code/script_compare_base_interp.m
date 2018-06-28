%% pick which data to interpolate here.
load('base_data.mat')

%% start the figure and plot the original data
figure;
subplot(1,2,1)
plotVField(xPos,yPos,repmat(2,size(xPos)),xVel,yVel,zVel,velMag);
xlim([4100,4800]); ylim([1.105e4,1.18e4]);
title('base scheme')



%% make the interpolated plot
subplot(1,2,2)
% pick interpolation resolution (half meters) and get velocity components
% on that grid
[X,Y] = meshgrid(4000:0.5:4800,1.1e4:0.5:1.18e4);
xVinterp = griddata(xPos,yPos,xVel,X,Y,'natural');
yVinterp = griddata(xPos,yPos,yVel,X,Y,'natural');
zVinterp = griddata(xPos,yPos,zVel,X,Y,'natural');
% compute the magnitude of velocity vectors to color and plot stuff
interpVals = sqrt(xVinterp.^2+yVinterp.^2+zVinterp.^2);
[XX,YY,ZZ] = meshgrid(1:size(X,1),1:size(Y,2),2);
plotVField(XX,YY,ZZ,xVinterp,yVinterp,zVinterp,interpVals);
title('interpolated base scheme')
xlim([100,size(X,1)]); ylim([50,size(Y,1)]);

%% plot the streamlines
figure; hold on;
% define the resolution
[startX,startY]=meshgrid(1:20:1601,1:20:size(Y,1));
% actually make them  (ignoring the z-velocity component)
streams = stream2(XX, YY, xVinterp,yVinterp,startX,startY);
streamline(streams);
xlim([100,size(X,1)]); ylim([50,size(Y,1)]);