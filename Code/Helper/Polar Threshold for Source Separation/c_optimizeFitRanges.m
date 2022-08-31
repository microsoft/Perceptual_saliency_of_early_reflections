clear all
close all

% Use this script to calculate the fitting interval with the minimal fit
% error for c_calculatePolarThreshold.m. 
% If the grid remains unchanged, this will yield the preset value of 114°.


%% Load results from baumgartner model:
pathname = './_Results/baumgartner';

filename = 'baumgartner_all_VPs';
checksum = 'a396ae7cf8afe4c00387625fa804ebe0'; % checksum of grid configuration
if ~exist(filename, 'var')
    fprintf('loading baumgartner results...')
    load(fullfile(pathname, [filename '_' checksum]));
    fprintf(' Done.\n')
end

%% Set parameters for polar threshold calculation:
threshold_W1 = 30.2332; 

RegAngles_W1 = 30:5:180;

%% Set other parameters:
% doPlot = 'polarThreshold_avg';
doPlot = false;
savePlot = sprintf('./_Plots (fit range optimization)/W1threshold%.4f', threshold_W1);
printProgress = true;

if savePlot
    if ~isfolder(savePlot)
        mkdir(savePlot)
    end
end

%% Allocate space 

polarThreshold_W1 = nan(numel(baumgartner_all_VPs.pol_ref), ...
                        numel(baumgartner_all_VPs.lat_ref), ...
                        2, numel(RegAngles_W1));
fitError_W1 = nan(numel(baumgartner_all_VPs.pol_ref), ...
                        numel(baumgartner_all_VPs.lat_ref), ...
                        2, numel(RegAngles_W1));

p = baumgartner_all_VPs.p;
pol = baumgartner_all_VPs.pol;
pol_delta = baumgartner_all_VPs.pol_delta;
pol_ref = baumgartner_all_VPs.pol_ref;
lat_ref = baumgartner_all_VPs.lat_ref;


%% parfor loop
parfor ll = 1:numel(baumgartner_all_VPs.lat_ref)
    [p_W1, Err_W1] = calculateThreshold(p(:,:,:,ll),                    ...
                                        pol,                            ...
                                        pol_delta,                      ...
                                        pol_ref,                        ...
                                        lat_ref(ll),                    ...
                                        threshold_W1,                   ...
                                        RegAngles_W1,                   ...   
                                        doPlot, savePlot, printProgress ... 
                                       );

    polarThreshold_W1(:,ll,:,:) = p_W1;
    fitError_W1(:,ll,:,:) = Err_W1;
end


%% fitError processing:
% average over polar reference angle:
fitError_W1_processed = mean(fitError_W1,1);
% average over positive and negative opening angles:
fitError_W1_processed = mean(fitError_W1_processed,3);
fitError_W1_processed = squeeze(fitError_W1_processed);

% weighted average over the lateral angle
weights = cos(deg2rad(lat_ref));
weightedmean = (weights * fitError_W1_processed) / sum(weights);

% fit quadratic curve
f = fit(RegAngles_W1(:),weightedmean(:),'poly2');
[~, RegAngle_W1] = min(f(1:180));

%Plot
AKf(20,8)
plot(RegAngles_W1,weightedmean)
hold on
plot(RegAngles_W1, f(RegAngles_W1))
legend({'data','quadratic fit'})
xline(RegAngle_W1,':','label', ...
    sprintf('optimal fit range: $[0^{\\circ}-%i^{\\circ}]$',RegAngle_W1), ...
    'interpreter','latex','HandleVisibility','off')
xlabel('$W_1$ fit range $[^{\circ}]$','interpreter','latex')
ylabel('Mean squared fit error','interpreter','latex')
title('$W_1$ fit range optimization ($w_i=\cos(lat_i)$ weighting for lateral average)','interpreter','latex')


%% Sanity check: Scalar value for reciprocity
% The polar threshold is not necessarily reciprocal, meaning that for
% example the polar threshold for lat=0, pol = -25, Delta pol > 0 is not
% necessarily equal to the polar threshold for lat = 0, pol= + 25,         
% Delta pol < 0
%
% Reciprocity is desired in the frontal and dorsal directions. Here, we
% look at the areas with a great circle distance < 45 degrees from the
% points (lat=0, pol=0) and (lat=0, pol=180).
%
% For each reference point with pol<0 (or pol<180 for dorsal direction),
% the polar threshold for positive openening angles is compared to the
% polar threshold for negative opening angles at the point mirrored to
% pol=0 (pol=180 for dorsal direction). If both polar threshold values are
% similar, the polar threshold is reciprocal.

% Find the points in the frontal direction:
lat_front = 0;
pol_front = 0;
lat_back  = 0;
pol_back  = 180;
maxdist = 45.1; %.1 is for rounding errors, otherwise <= maxdist will get unexpected results

%distance.m uses lateral and longitudinal angles (geographic coordinate
%system) whereas we use the lateral-polar angle coordinate system. if we
%switch the lat and pol angle we get the same coordinate system as the
%geographig coordinate system, only rotated by 90 degrees. Thus, the grat
%circle distance will be calculated correctly.

reciprocity = nan(numel(baumgartner_all_VPs.pol_ref), ...
                        numel(baumgartner_all_VPs.lat_ref), ...
                        1, numel(RegAngles_W1));

% calculate indices of front and back polar angle:
idx_front = find(pol_ref == pol_front);
idx_back = find(pol_ref == pol_back);

for ll = 1:numel(lat_ref)
    for pp = 1:numel(pol_ref)
        % check if point is in the frontal direction 
        if pol_ref(pp) >= -maxdist && pol_ref(pp) < 0 && ...
            distance(pol_ref(pp),lat_ref(ll),pol_front,lat_front) <= maxdist
            %find index of the point mirrored at pol_front:
            idx = pp;
            idx_mirrored = idx_front + (idx_front - idx);
            
            % get polar thresholds and calculate ratio
            for rr = 1:numel(RegAngles_W1)
                % get polar threshold at reference point (pos opening angle):
                x = polarThreshold_W1(idx         ,ll,1,rr);
                % get polar threshold at mirrored point (neg opening angle):
                y = polarThreshold_W1(idx_mirrored,ll,2,rr);
                % calculate ratio of absolute values:
                reciprocity(pp,ll,1,rr) = max(abs(x),abs(y)) / min(abs(x),abs(y));
            end
            
        % check if point is in the dorsal direction 
        elseif pol_ref(pp) >= pol_back-maxdist && pol_ref(pp) < 180 && ...
            distance(pol_ref(pp),lat_ref(ll),pol_back,lat_back) <= maxdist
            %find index of the point mirrored at pol_back:
            idx = pp;
            idx_back = find(pol_ref == pol_back);
            idx_mirrored = idx_back + (idx_back - idx);
            
            % get polar thresholds and calculate ratio
            for rr = 1:numel(RegAngles_W1)
                % get polar threshold at reference point (pos opening angle):
                x = polarThreshold_W1(idx         ,ll,1,rr);
                % get polar threshold at mirrored point (neg opening angle):
                y = polarThreshold_W1(idx_mirrored,ll,2,rr);
                % calculate ratio of absolute values:
                reciprocity(pp,ll,1,rr) = max(abs(x),abs(y)) / min(abs(x),abs(y));
            end
        end
        
        
    end
end
%
% Average across polar angles
reciprocity = squeeze(mean(reciprocity,1,'omitnan'));

reciprocity(isnan(reciprocity)) = 0;
% Weighted average across lateral angles: 
reciprocity = (weights * reciprocity) / sum(weights);

% fit reciprocity:
f = fit(RegAngles_W1(:),reciprocity(:),'poly5');
reciprocity_opti = interp1(f(1:180),1:180,1);

%Plot
AKf(20,8)
plot(RegAngles_W1,reciprocity)
hold on
plot(30:180,f(30:180),'--')
legend({'data','quintic fit'})
yline(1,'--','HandleVisibility','off')
xline(reciprocity_opti,':','label','optimal reciprocity', ...
    'interpreter','latex','HandleVisibility','off')
xline(RegAngle_W1,':','label','optimal fit error', ...
    'interpreter','latex','HandleVisibility','off')

xlabel('$W_1$ fit range $[^{\circ}]$','interpreter','latex')
ylabel('Reciprocity','interpreter','latex')
title('$W_1$ reciprocity ($w_i=\cos(lat_i)$ weighting for lateral average)','interpreter','latex')



