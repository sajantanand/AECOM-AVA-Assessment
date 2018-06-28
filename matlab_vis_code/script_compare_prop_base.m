% Plot just the vecotr fields side-by-side with vectors colored by
% magnitude or TKE
figure; hold on;

load('proposal_data.mat')
% subplot(1,2,1)
plotVField(xPos,yPos,repmat(2,size(xPos)),xVel,yVel,zVel,velMag);
title('proposed scheme')

% % load('FULL_proposal_data.mat')
% % % % subplot(1,2,1)
% % cond = (xcoordinate>3601 & xcoordinate <5000);
% % x = xcoordinate(cond);
% % y = ycoordinate(cond);
% % z= zcoordinate(cond);
% % u = xvelocity(cond);
% % v = yvelocity(cond);
% % w = zvelocity(cond);
% % mags = velocitymagnitude(cond);
% % cond = (y>1.05e4 & y<1.25e4);
% % x = x(cond);
% % y = y(cond);
% % z = z(cond);
% % u = u(cond);
% % v = v(cond);
% % w = w(cond);
% % mags = mags(cond);
% % plotVField(x,y,z,u,v,w,mags);
% % xlim([4000,4800])
% % ylim([1.1e4,1.18e4])
% % zlim([-5,200]);
% % title('proposed scheme')

load('base_data.mat')
subplot(1,2,2)
plotVField(xPos,yPos,repmat(2,size(xPos)),xVel,yVel,zVel,velMag);
title('base scheme')

