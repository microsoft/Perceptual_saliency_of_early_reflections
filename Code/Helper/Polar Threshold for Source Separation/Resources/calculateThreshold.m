function [polarThreshold_W1, fitError_W1] =                             ...
            calculateThreshold(p, pol, pol_delta,                       ...
                               pol_ref, lat_ref,                        ...
                               threshold_W1,                            ...
                               RegAngles_W1,                            ...
                               doPlot, savePlot, printProgress)
%calculateThreshold Calculates polar Threshold from baumgartner model
% Written by Tobias Jüterbock
% See calculateBaumgartner.m for more detailed information on what happens
% there. p is a 4-D array with probability distributions of the perceived
% incident angles for summed HRTFs of two angles of incidents. 
%
% The baumgartner model is designed for polar localization estimates of one
% incident HRTF and not suitable for analyzing two simultaneous events.
% Nevertheless, the baumgartner model yields reasonable results if the
% events are spatially close to each other (usually two sources close to
% each other are interpreted as one phantom source located somewhere
% between the actual source positions). If the sources are not close to
% each other, the model does not yield reasonable results anymore.
%
% The idea of this threshold is finding the minimal opening angle where the
% perception of two sources is distinguishable from one source being there.
% We assume that this is the case when the baumgartner model starts falling
% apart and stops interpreting two summed sources as one source. 
%
% The probability distributions are analyzed and fitted with gauss
% distributions. The goodness of the fit (Rsquare) is then used as an
% indicator for the model interpreting the input as one source. The output
% is the first opening angle where Rsquare falls below the value specified
% in 'threshold'.
%
% I N P U T:
% p        : probability distributions from the baumgartner model as a 4-D
%            array of size P x D x R x L
% pol      : polar angles to be tested (size: P x 1)
% pol_delta: polar opening angles (size: 1 x D)
% pol_ref  : polar reference angles (size: 1 x R)
% lat_ref  : lateral reference angles (size: 1 x L)
% threshold: threshold for goodness of fit


% Reminder: dimensions of p
% 1st dim: polar angles to be tested
% 2nd dim: delta_polar angles to be tested
% 3rd dim: polar reference angle
% 4th dim: lateral reference angle
      
%% Parameters:

%% Check directory for Plots
plotTitle = '';
if ischar(doPlot)
    plotTitle = doPlot;
    doPlot    = true;
end

plotpath = '.';
if ischar(savePlot)
    plotpath = savePlot;
    savePlot = true;
end

% create Plot directories if necessary:
if ~isfolder(fullfile(plotpath))
    mkdir(fullfile(plotpath));
end
if ~isfolder(fullfile(plotpath, 'fig'))
    mkdir(fullfile(plotpath,'fig'));
end
if ~isfolder(fullfile(plotpath,'pdf')) 
    mkdir(fullfile(plotpath,'pdf'));
end
if ~isfolder(fullfile(plotpath,'eps')) 
    mkdir(fullfile(plotpath,'eps'));
end


%% allocate space for results:
% 1st dim: polar reference angle
% 2nd dim: lateral reference angle
% 3rd dim: 1 for positive opening angles, 2 for negative opening angles
polarThreshold_W1    = nan(numel(pol_ref), numel(lat_ref), 2, numel(RegAngles_W1));
fitError_W1          = nan(numel(pol_ref), numel(lat_ref), 2, numel(RegAngles_W1));
fitparams_W1         = cell(numel(pol_ref), numel(lat_ref), 2, numel(RegAngles_W1));

%% Loop over all test conditions:
%% Loop over lateral reference angles:
for ll = 1:numel(lat_ref)
    if printProgress
       fprintf('Calculating lat %i of %i\n',ll,numel(lat_ref)) 
    end
    % Calculate threshold for current lateral angle
    threshold_W1_lat = thresholdLateralIncrease(threshold_W1,lat_ref(ll));
    
    %% Loop over polar reference angles:
    for pp = 1:numel(pol_ref)
        if printProgress
            fprintf('Calculating lat %i of %i, pol %i of %i\n', ...
                    ll,numel(lat_ref), pp,numel(pol_ref)) 
        end
        % reset W1 values:
        W1 = nan(size(pol_delta));
        
        PDFs = p(:,:,pp,ll);
        % Reslice the data so the values center around pol_ref:
        % find index of pol at pol_ref+180:
        id = find (pol == mod(pol_ref(pp)+180+90, 360)-90);
        if id ~= 360
            PDFs = [PDFs(id+1:end,:);PDFs(1:id,:)];
            pol2 = (pol(id+1):pol(id)+360).';            
        else
            pol2 = pol;
        end
        
        % Convert polar reference angle to the resliced polar angles:
        pol_ref2 = pol_ref(pp);
        while ~(pol_ref2 >= pol2(1) && pol_ref2 <= pol2(end))
            pol_ref2 = pol_ref2 + 360;
        end
        
        %% Loop over polar opening angles:
        for dd = 1:numel(pol_delta)
            %% calculate the Wasserstein Distance:
            while isnan(W1(dd))
                % sample random numbers following the PDFs
                h1 = randpdf(PDFs(:,1 ), pol2, [10000,1]);
                h2 = randpdf(PDFs(:,dd), pol2, [10000,1]);
                % calculate the Wasserstein Distance from the samples:
                W1(dd) = ws_distance(h1,h2,1);
                % ws_distance can return NaN sometimes, in those cases the
                % while loop is necessary.
            end
        end
        
        %% Calculate the polar Threshold for positive opening angles:
        % separate positive (0 to 180) and negative (180 to 360) opening
        % angles,
        % first, we only consider positive opening angles:
        
        %% Calculate the polar threshold with regression fit to W1:
        for i = 1:numel(RegAngles_W1)
            % Perform linear fit and calculate intersection with threshold:
            [polarThreshold_W1(pp,ll,1,i),                              ...
             fitparams_W1{pp,ll,1,i}                                    ...
            ] = linearFitAndThreshold(pol_delta,W1,                     ...
                                      threshold_W1_lat,                 ...
                                      [0, RegAngles_W1(i)],             ...
                                      180);
        end

        %% Calculate the polar Threshold for negative opening angles:
        % Calculate the polar Threshold with regression fit to W1:
        for i = 1:numel(RegAngles_W1)
            % Perform linear fit and calculate intersection with threshold:
            [polarThreshold_W1(pp,ll,2,i),                              ...
             fitparams_W1{pp,ll,2,i}                                    ...
            ] = linearFitAndThreshold(-pol_delta, flip(W1),             ...
                                      threshold_W1_lat,                 ...
                                      [0, -RegAngles_W1(i)],            ...
                                      -180);
        end
   
        %% Plot baumgartner model and thresholds
        if doPlot
            plotcolors = [[0.4660, 0.6740, 0.1880]; ...
                          [0.9290, 0.6940, 0.1250]; ...
                          [0.8500, 0.3250, 0.0980];...
                          [0.6350, 0.0780, 0.1840];	
                          [0.4940, 0.1840, 0.5560]];...
                          
            %% Plot Results for each reference point
            AKf(12,7.5)
            t = tiledlayout(2,1);
            t.Padding = 'compact';
            nexttile
            imagesc([pol_delta,pol_delta(1)+360] , pol2, [PDFs , PDFs(:,1)])
            hold on
            plot(pol_delta, mod(pol_ref(pp)-pol2(1),360)+pol2(1).*ones(size(pol_delta)), 'color', [AKcolors('y') .5], 'LineWidth', 1.5)
            plot(pol_delta, mod(pol_ref(pp).*ones(size(pol_delta))+pol_delta-pol2(1), 360)+pol2(1), 'color', [AKcolors('y') .5], 'LineWidth', 1.5)
            colorbar
            xlim([0, 360])
            ylim([pol2(1),pol2(end)])
            xticks(0:90:360)
            tick = unique(sort(mod([-90, 0, 90, 180, 270]+360-pol2(1),360)+pol2(1)));
            yticks(tick);
            yticklabels(mod(tick+90,360)-90);
            set(gca,'yDir','normal')
            if ~isnan(polarThreshold_W1(pp,ll,1,:))
                xline(polarThreshold_W1(pp,ll,1,1), 'w:');
                if size(polarThreshold_W1,4) >= 2
                xline(polarThreshold_W1(pp,ll,1,2), 'w:');
                end
                if size(polarThreshold_W1,4) >= 3
                xline(polarThreshold_W1(pp,ll,1,3), 'w:');
                end
            end
            if ~isnan(polarThreshold_W1(pp,ll,2,:))
                xline(360+polarThreshold_W1(pp,ll,2,1),'w:');
                if size(polarThreshold_W1,4) >= 2
                xline(360+polarThreshold_W1(pp,ll,2,2),'w:');
                end
                if size(polarThreshold_W1,4) >= 3
                xline(360+polarThreshold_W1(pp,ll,2,3),'w:');
                end
            end
            
            ylabel('$\theta$ [$^{\circ}$]','interpreter','latex')
          
            nexttile
            plot([pol_delta,pol_delta(1)+360], [W1, W1(1)],'HandleVisibility','off');
            hold on
            for i = 1:numel(RegAngles_W1)
                % plot W1 regression:
                W1RegAngle = RegAngles_W1(i);
                xmax = max(W1RegAngle, polarThreshold_W1(pp,ll,1,i));
                semilogy(0:xmax, (fitparams_W1{pp,ll,1,i}(0:xmax)),'color',plotcolors(i,:));
                xline(polarThreshold_W1(pp,ll,1,i),':','HandleVisibility','off')
            end
            for i = 1:numel(RegAngles_W1)
                % plot W1 regression:
                W1RegAngle = -RegAngles_W1(i);
                xline(360+polarThreshold_W1(pp,ll,2,i),':','HandleVisibility','off')
                xmin = min(W1RegAngle, polarThreshold_W1(pp,ll,2,i));
                semilogy(360+(xmin:0), (fitparams_W1{pp,ll,2,i}(xmin:0)),'color',plotcolors(i,:))
            end
            xlim([0 360])
            ylim([0 120])
            xticks(0:90:360)
            yl = yline(threshold_W1_lat, '--','$t_{W_1}$','labelVerticalAlignment','top','interpreter','latex');
            yl.LabelHorizontalAlignment = 'left';
            if numel(RegAngles_W1) > 1
            legend('0-'+string(RegAngles_W1)+'°');
            legend('boxoff')
            else 
                legend('linear fit')
            end
            xlabel('$\Delta \theta$ [$^{\circ}$]','interpreter','latex')
            ylabel('$W_1$', 'interpreter','latex')
            hold off
            
            if savePlot
                filename = sprintf('%s_lat%i_pol%i.fig',plotTitle,lat_ref(ll),pol_ref(pp));
                savefig(fullfile(plotpath,'fig',filename))
                filename = sprintf('%s_lat%i_pol%i.pdf',plotTitle,lat_ref(ll),pol_ref(pp));
                saveas(gcf, fullfile(plotpath,'pdf',filename))
                filename = sprintf('%s_lat%i_pol%i.eps',plotTitle,lat_ref(ll),pol_ref(pp));
                saveas(gcf, fullfile(plotpath,'eps',filename))
                close(gcf)
            end


        end
        
        %% Calculate fit error:
        % Calculate fit error for W1:
        for i = 1:numel(RegAngles_W1)
            % Positive opening angles:
            % get W1 values from 0 to polar threshold
            x = W1(pol_delta <= polarThreshold_W1(pp,ll,1,i));
            x = x(:);
            % get fit values from 0 to polar threshold
            y = fitparams_W1{pp,ll,1,i};
            y = y(pol_delta(pol_delta <= polarThreshold_W1(pp,ll,1,i)));
            y = y(:);
            % calculate mean squared error:
            fitError_W1(pp,ll,1,i) = immse(x,y);
            
            % Negative opening angles:
            % get W1 values from 0 to polar threshold
            x = W1(pol_delta >= polarThreshold_W1(pp,ll,2,i)+360);
            x = x(:);
            % get fit values from 0 to polar threshold
            y = fitparams_W1{pp,ll,2,i};
            y = y(pol_delta(pol_delta >= polarThreshold_W1(pp,ll,2,i)+360)-360);
            y = y(:);
            % calculate mean squared error:
            fitError_W1(pp,ll,2,i) = immse(x,y);
        end
    end
end
end