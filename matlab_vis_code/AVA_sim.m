%% Current AVA practices
% In the current practice, we test points along the perimeter and points
% around the building in a radius equal to the height of the tallest
% building. Here we demonstrate this. 
% Set show_points to 0 to not display the points picked or 1 to see them
%% Define paramaters/ initialization
% Display test points whenever base scheme is better than proposed?
show_bad_pts = 0;

% Need to initialize or already done?
if ~exist('dont_have_data_and_masks','var')
  dont_have_data_and_masks = 1;
else
  dont_have_data_and_masks = 0;
end

% Params for results and testing
test_pt_iterations = 100;
% 30-50 perimeter points usually
default_num_perimeter_points = 10;
% 50-80 total points usually
default_num_total_points = 10;
% space the points about 10-50 meters apart
ptSpacing = 20; %10m apart
bins = 20;  % for histograms comparing velocity ratios
plot_results = 1;

% Show each time test points are picked? (should really only be on with
% very few iterations to make a figure for a presentation or to see what
% the configurations look like
show_points = 0; % always shows first set of points though
% probably shouldn't remove the next line.  It catches when you accidently
% set iterations very high but left show_points on.
show_points = show_points && (test_pt_iterations < 4);

% setup empty variables
p_SVR = zeros(test_pt_iterations,1);
p_LVR = zeros(test_pt_iterations,1);
b_SVR = zeros(test_pt_iterations,1);
b_LVR = zeros(test_pt_iterations,1);
p_bad_SVR = [];
b_bad_SVR = [];
p_bad_LVR = [];
b_bad_LVR = [];
%% Data setup and preprocessing
if dont_have_data_and_masks 
  
  %% %% Step 1: Retrieve the data and interpolate it to a uniform gridspace 
  %% Step 1a: Interpolate base scheme
  load('full_site_data.mat');

  [X,Y] = meshgrid(4000:0.5:4800,1.1e4:0.5:1.18e4);
  xVinterp = griddata(xPos,yPos,xVel,X,Y,'natural');
  yVinterp = griddata(xPos,yPos,yVel,X,Y,'natural');
  zVinterp = griddata(xPos,yPos,zVel,X,Y,'natural');
  interpVals_base1 = sqrt(xVinterp.^2+yVinterp.^2+zVinterp.^2);
  b_interp = zeros(size(X));
  for ii = 1:size(X,1)
    b_interp(size(X,1)-ii+1,:) = interpVals_base1(ii,:);
  end
  clear('interpVals_base1');

  %% Step 1b: Interpolate proposal scheme
  load('proposal_data.mat');

  [X,Y] = meshgrid(4000:0.5:4800,1.1e4:0.5:1.18e4);
  xVinterp = griddata(xPos,yPos,xVel,X,Y,'natural');
  yVinterp = griddata(xPos,yPos,yVel,X,Y,'natural');
  zVinterp = griddata(xPos,yPos,zVel,X,Y,'natural');
  interpVals_prop1 = sqrt(xVinterp.^2+yVinterp.^2+zVinterp.^2);
  p_interp = zeros(size(X));
  for ii = 1:size(X,1)
    p_interp(size(X,1)-ii+1,:) = interpVals_prop1(ii,:);
  end
  clear('interpVals_prop1','xVinterp','yVinterp','zVinterp','xPos','yPos','xVel','yVel','zVel','X','Y');

  %% %% Step 2: Specify regions for sampling
  % Select a polygon hugging the development schemes
  % Points are to be tested within a radius of H around the development, so
  % select this region.  It's also necessary to specify the height of the
  % buildings.  Here since the *scale is .5 m per unit* in the interpolation,
  % we multiply the height input by 2.

  %% Step 2a: Scheme masks
  per_masked_base = b_interp;
  building_mask = b_interp <.01 | p_interp <.01;
  
  roi_mask =roipoly(building_mask);
  close(gcf)
  
  % get mask for building perimeter (SVR)
  perim_rad = 10;
  buildings = (b_interp <.01) & roi_mask;
  % per_mask is just around the perimeter of the development and not in the development
  per_mask = (~building_mask) & (~buildings) & (conv2(double(roi_mask),double(fspecial('disk',2*perim_rad)>0),'same')>0);
  per_mask =  per_mask & ~roi_mask; 
  % get mask for radius 98.5m area around the proposed region (for LVR)
  height = 98.5;%(in meters)
  h_mask = (~per_mask) & (~building_mask) & conv2(double(roi_mask),double(fspecial('disk',2*height)>0),'same')>0;
  figure; hold on; imshow(h_mask); title('base H mask'); hold off
  figure; hold on; imshow(per_mask); title('perimeter mask'); hold off
end
%% Test Points
progress_bar = waitbar(0,'Please wait...');
for iter = 1:test_pt_iterations
  %% %% Step 3: Trials to compute SVR,LVR
  % print progress every 10% of total iterations
 if mod(iter,round(.01*test_pt_iterations))==0
    frac = iter/test_pt_iterations;
    percent = round(frac*100);
    waitbar(frac, progress_bar, [num2str(percent),'% total trials completed...']);
 end
  %% %% Step 3 a1: Pick perimeter points proposed scheme and calculate SVR, LVR
  % 30 - 50 perimeter points
  % 50 - 80 overall test points
  % space the points about 10-50 meters apart
  
  numPerPts = default_num_perimeter_points;
  numTotPts = default_num_total_points;
  per_pts = zeros(numPerPts,2);
  remPerMask = per_mask;
  inds = zeros(numTotPts,1);
  % select the perimeter points
  for ii = 1:numPerPts
    index  = find(remPerMask);
    if ~isempty(index)
      select = index(randperm(length(index), 1));
      if p_interp(select)>0 && b_interp(select)>0
        inds(ii) = select;
        % next line really should be a circle, but it removes a region from
        % possible candidates (enforces point spacing constraint)
        [y,x] = ind2sub(size(remPerMask), select);
        per_pts(ii,:) = [y,x];
        remPerMask(y-ptSpacing:y+ptSpacing,x-ptSpacing:x+ptSpacing) = zeros(2*ptSpacing+1,2*ptSpacing+1);
        figure; imshow(remPerMask)
      else
        ii = ii-1;
      end
    else
      numPerPts = ii-1;
      per_pts = per_pts(1:numPerPts,:);
      break;
    end
  end

  %% %% Step 3 a2: Pick non-perimeter points proposed scheme and calculate SVR, LVR
  numNonPer = numTotPts-numPerPts;
  Pts = zeros(numNonPer,2);
  remNonPerMask = (h_mask & ~per_mask);
  % select non-perimeter points
  for ii = 1:numNonPer
    index  = find(remNonPerMask);
    if ~isempty(index)
      select = index(randperm(length(index), 1));
      if p_interp(select)>0 && b_interp(select)>0
        inds(ii+numPerPts) = select;
        % next line really should be a circle, but it removes a region from
        % possible candidates (enforces point spacing constraint)
        [y,x] = ind2sub(size(remNonPerMask), select);
        Pts(ii,:) = [y,x];
        remNonPerMask(y-ptSpacing:y+ptSpacing,x-ptSpacing:x+ptSpacing) = zeros(2*ptSpacing+1,2*ptSpacing+1);
      else
        ii = ii-1;
      end
    else
      numNonPer = ii-1;
      Pts = Pts(1:numNonPer,:);
      break;
    end
  end
  numTotPts = numNonPer + numPerPts;
  inds = inds(1:numTotPts);
  %% %% Step 4: Calculate SVR, LVR Show the proposed points
  % SVR is average of all the perimeter points
  % LVR is average of all of the points
  p_per_vals = p_interp((inds(1:numPerPts)));
  p_all_vals = p_interp((inds));
  
  b_per_vals = b_interp((inds(1:numPerPts)));
  b_all_vals = b_interp((inds));
  
  p_SVR(iter) = mean(p_per_vals);
  p_LVR(iter) = mean(p_all_vals);
  
  b_SVR(iter) = mean(b_per_vals);
  b_LVR(iter) = mean(b_all_vals);
  if show_points || (iter == 1)
    plot_test_pts(p_interp,per_pts,Pts,height);
    title(['Proposed scheme with testpoints, SVR: ',num2str(p_SVR(iter)),', LVR: ',num2str(p_LVR(iter))]);
    hold off
    plot_test_pts(b_interp,per_pts,Pts,height);
    title(['Base scheme with testpoints, SVR: ',num2str(b_SVR(iter)),', LVR: ',num2str(b_LVR(iter))]);
    hold off
  end
  
  %% %% 5 Check out bad configs
  if (b_LVR(iter)>p_LVR(iter) || b_SVR(iter)>p_SVR(iter))
    if show_bad_pts
      warning('--proposed scheme failed test--')

      % base test points
      plot_test_pts(b_interp,b_perPts,b_Pts,height)
      title(['(BAD-config) Base scheme with testpoints, SVR: ',num2str(b_SVR(iter)),', LVR: ',num2str(b_LVR(iter))]);
      hold off

      figure; hold on
      histogram(b_all_vals,bins);
      line([b_LVR(iter), b_LVR(iter)], ylim, 'LineWidth', 2, 'Color', 'b');
      histogram(b_per_vals,bins);
      line([b_SVR(iter), b_SVR(iter)], ylim, 'LineWidth', 2, 'Color', 'r');
      legend('all points','LVR','perimeter points','SVR')
      title('speed at all test points in base scheme')
      hold off

      % prop test points
      plot_test_pts(p_interp,per_pts,Pts,height);
      title(['(BAD-config) Proposed scheme with testpoints, SVR: ',num2str(p_SVR(iter)),', LVR: ',num2str(p_LVR(iter))]);
      hold off

      figure; hold on
      histogram(p_all_vals,bins*2);
      line([p_LVR(iter), p_LVR(iter)], ylim, 'LineWidth', 2, 'Color', 'b');
      histogram(p_per_vals,bins);
      line([p_SVR(iter), p_SVR(iter)], ylim, 'LineWidth', 2, 'Color', 'r');
      legend('all points','LVR','perimeter points','SVR')
      title('speed at all test points in proposed scheme')
      hold off
    end   
    % update the list of failed cases
    b_bad_SVR = [b_bad_SVR;b_SVR(iter)];
    b_bad_LVR = [b_bad_LVR;b_LVR(iter)];
  end
% %     ['num too high pts prop (LVR +2std): ' ,num2str(sum(p_pt_vals>mean(p_pt_vals)+2*sqrt(var(p_pt_vals))))]
% %     ['num too low pts prop (LVR -2std): ' ,num2str(sum(p_pt_vals>mean(p_pt_vals)+2*sqrt(var(p_pt_vals))))]
  p_bad_SVR = [p_bad_SVR;p_SVR(iter)];
  p_bad_LVR = [p_bad_LVR;p_LVR(iter)];
end
close(progress_bar);
%% Show the results
mlp = mean(p_LVR);
slp = sqrt(var(p_LVR));
mlb = mean(b_LVR);
slb = sqrt(var(b_LVR));

msp = mean(p_SVR);
ssp = sqrt(var(p_SVR));
msb = mean(b_SVR);
ssb = sqrt(var(b_SVR));

if plot_results
  %% Step 6: plot histograms + gaussian fitof LVR and SVR comparisons
  figure; hold on
  hsb =histfit(b_SVR,bins,'normal'); hsb(1).FaceAlpha = .6; hsb(1).FaceColor = [.1 .6 .8]; hsb(2).Color = 'b';
  hsp =histfit(p_SVR,bins,'normal'); hsp(1).FaceAlpha = .6; hsp(1).FaceColor = [.8 .6 .1]; hsp(2).Color = 'r';
  title(['SVR over ',num2str(test_pt_iterations),' trials of random test points'])
  legend('base SVR','gfit base','prop SVR','gfit prop')
  hold off

  figure; hold on
  hlb =histfit(b_LVR,bins,'normal'); hlb(1).FaceAlpha = .6; hlb(1).FaceColor = [.1 .6 .8]; hlb(2).Color = 'b';
  hlp =histfit(p_LVR,bins,'normal'); hlp(1).FaceAlpha = .6; hlp(1).FaceColor = [.8 .6 .1]; hlp(2).Color = 'r';
  title(['LVR over ',num2str(test_pt_iterations),' trials of random test points'])
  legend('base LVR','gaussian fit base','prop LVR','gaussian fit prop')
  hold off

  figure; hold on
  hs =histfit(p_SVR-b_SVR,bins,'normal'); hs(1).FaceAlpha = .6; hs(1).FaceColor = [.1 .6 .8]; hs(2).Color = 'b';
  hl =histfit(p_LVR-b_LVR,bins,'normal'); hl(1).FaceAlpha = .6; hl(1).FaceColor = [.8 .6 .1]; hl(2).Color = 'r';
  title(['differences in SVR and LVR over ',num2str(test_pt_iterations),' trials of random test points'])
  legend('prop minus base (SVR)','gaussian fit SVR diff','prop minus base (LVR)','gaussian fit LVR diff')
  hold off
end
