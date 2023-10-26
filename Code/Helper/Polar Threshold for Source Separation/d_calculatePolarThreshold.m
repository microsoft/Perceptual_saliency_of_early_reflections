close all; clear; clc
% requires Parallel Computing Toolbox. If not available, change the parfor
% loop into a for loop


%% Load results from baumgartner model:
pathname = './_Results/baumgartner';
filename = 'baumgartner_all_VPs';
checksum = 'a396ae7cf8afe4c00387625fa804ebe0'; % checksum of grid configuration


if ~exist(filename, 'var')
    fprintf('loading baumgartner results...')
    load(fullfile(pathname, [filename,'_',checksum]));
    fprintf(' Done.\n')
end

%% Set parameters for polar threshold calculation:
threshold_W1 = 30.2332; 
RegAngles_W1 = 114;

%% Set other parameters:
doPlot = 'polarThreshold_avg';
savePlot = sprintf('./Plots (full grid)/W1threshold%.4f', threshold_W1);
printProgress = true;
saveResults = true;

%%
polarThreshold_W1 = nan(numel(baumgartner_all_VPs.pol_ref), ...
                        numel(baumgartner_all_VPs.lat_ref), ...
                        2, numel(RegAngles_W1));
p         = baumgartner_all_VPs.p;
pol       = baumgartner_all_VPs.pol;
pol_delta = baumgartner_all_VPs.pol_delta;
pol_ref   = baumgartner_all_VPs.pol_ref;
lat_ref   = baumgartner_all_VPs.lat_ref;
parfor ll = 1:numel(baumgartner_all_VPs.lat_ref)
    p_W1 = calculateThreshold(p(:,:,:,ll),                              ...
                              pol,                                      ...
                              pol_delta,                                ...
                              pol_ref,                                  ...
                              lat_ref(ll),                              ...
                              threshold_W1,                             ...
                              RegAngles_W1,                             ...   
                              doPlot, savePlot, printProgress           ... 
                             );

    polarThreshold_W1(:,ll,:,:) = p_W1;
end

%% Save results
if saveResults == true
    filename_results = ['polarThreshold_all_VPs','_',checksum,'.mat'];
    pathname_results = './Results/Threshold';
    if ~isfolder(pathname_results)
        mkdir(pathname_results)
    end
    polarThreshold = polarThreshold_W1(:,:,:,2);
    save(fullfile(pathname_results, filename_results), ...
         'lat_ref', 'pol_ref', 'polarThreshold');
end

%% Plot resullts
AKf(10*numel(RegAngles_W1),20)
idx = 1;

for jj = 1:2
    for ii  = 1:numel(RegAngles_W1)        
        subplot(2, numel(RegAngles_W1), idx)
        imagesc(baumgartner_all_VPs.lat_ref, ...
        baumgartner_all_VPs.pol, ...
        polarThreshold_W1(:,:,jj,ii))
        set(gca,'yDir','normal')
        
        if jj == 1
        AKpSetColormapAndColorbar_gca('Reds', 'EastOutside', 20, [0 180], 'polar Threshold')
        title(sprintf('fit range 0-%i°, pos. opening angle', RegAngles_W1(ii)))
        elseif jj == 2
        AKpSetColormapAndColorbar_gca('Reds_flip', 'EastOutside', 20, [-180 0], 'polar Threshold')
        title(sprintf('fit range: 0-%i°, neg. opening angle', RegAngles_W1(ii)))
        end

        xticks(-90:30:90)
        yticks(-90:90:270)
        idx = idx + 1;
    end
end
for i = [1,numel(RegAngles_W1)+1]
    subplot(2, numel(RegAngles_W1), i)
    ylabel('pol_{ref}')
end
for i = (numel(RegAngles_W1)+1):(numel(RegAngles_W1)*2)
    subplot(2, numel(RegAngles_W1), i)
    xlabel('lat_{ref}')
end
sgtitle(sprintf('polar threshold, linear regression of W_1, threshold: %.1f', threshold_W1))

if ischar(savePlot)
   plotpath = savePlot;
   savePlot = true;
end

if savePlot
    filename = sprintf('polarThreshold_W1');
    savefig(fullfile(plotpath,filename))
    saveas(gcf, fullfile(plotpath,[filename,'.pdf']))
end