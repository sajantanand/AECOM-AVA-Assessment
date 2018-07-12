function [] = plot_test_pts(wind_speed,perPts,otherPts,H)
  figure; imagesc(wind_speed); hold on
  scatter(perPts(:,2), perPts(:,1), 'mo','filled','MarkerEdgeColor','white', 'SizeData', 17);
  scatter(otherPts(:,2),otherPts(:,1),'ro','filled','MarkerEdgeColor',[.7,0,0],'SizeData',17);
  cc = parula; cc(1,:)=[0,0,0];
  legend([num2str(size(perPts,1)),' perimeter points'],[num2str(size(otherPts,1)),' points within H=',num2str(H),' meters']); 
  colormap(cc);
  c = colorbar;
  c.Label.String = 'Wind speed (m/s)';
  xlim([200,1400]); ylim([200,1400]);
end