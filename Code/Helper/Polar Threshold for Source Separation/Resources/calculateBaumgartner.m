function data = calculateBaumgartner(lat_ref, pol_ref, pol_mu, VP, printProgress, doPlot)
%CALCULATEBAUMGARTNER Calculates localization based on baumgartner model
% needs auditory modeling toolbox and AKHrirDatabase/HUTUBS Database
% Written by Tobias Jüterbock
% 
% I N P U T:
% lat_ref      : lateral reference angles (size: 1 x L)
% pol_ref      : polar reference angles (size: 1 x P)
% pol_mu       : stepsize in degrees for opening angle between sources
% VP           : Test subject whose HRTFs will be used (see AKHrirDatabase)
% printProgress: true to output progress to console
% doPlot       : true to plot results

% O U T P U T:
% data : struct with all relevant fields:
% data.p : 4D-array with the probability distributions from the baumgartner
%          model (size: 360 x 360/pol_mu x P x L)
% data.pol      : polar angles (size: 1 x 360);
% data.pol_mu   : pol_mu;
% data.pol_delta: opening angles (size: 1 x 360/pol_mu);

%% Input parameter check:
assert(sign(pol_mu) == 1, 'pol_mu has to be positive');

%% Allocate space 
VPstring = sprintf('pp%i_simulated_sh',VP);

pol = (-90:1:269)'; % polar angles at which reference HRTFs are provided
pol_delta = 0:abs(pol_mu):359;

data = struct;
data.pol = pol;
data.pol_mu = pol_mu;
data.pol_delta = pol_delta;
data.lat_ref = lat_ref;
data.pol_ref = pol_ref;
data.model = 'baumgartner';

% Allocate space for results:
% 1st dim: polar angles to be tested
% 2nd dim: delta_polar angles to be tested
% 3rd dim: polar reference angle
% 4th dim: lateral angle
baumgartner = nan(numel(pol), numel(pol_delta), numel(pol_ref), numel(lat_ref), 'single');

%% get the common transfer function
sg = sofia_fliege(256, false);
sg = [sg(:,1) pi/2-sg(:,2)] /pi*180;

H = AKhrirDatabase(VPstring, sg, false, 'SOFA');
H.GLOBAL_Comment = '';

[~, CTF] = SOFAhrtf2dtf(H);

clear H sg

%% loop across test conditions:
for ll = 1:numel(lat_ref)
    if printProgress
       fprintf('Calculating lat %i of %i\n',ll,numel(lat_ref)) 
    end
    for pp = 1:numel(pol_ref)
        %% get HRIRs
        % reference directional HRIRs
        POL_ref       = (-90:1:269)'; % polar angles at which reference HRTFs are provided
        [az, el]      = hor2sph(lat_ref(ll) * ones(size(POL_ref)), POL_ref);
        H_ref         = AKhrirDatabase(VPstring, [az el], false, 'SOFA');
        H_ref.Data.IR = ifft(fft(H_ref.Data.IR, [], 3) ./ fft(CTF.Data.IR, [], 3), [], 3, 'symmetric');

        clear az el

        % summed directional HRIRs
        POL_test      = (pol_ref(pp)+pol_delta)';     % polar angles for testing
        POL_test      = mod(POL_test+90, 360)-90; % map to -90 <= pol <= 270
        [az, el]      = hor2sph(lat_ref(ll) * ones(size(POL_test)), POL_test);
        H_sum         = AKhrirDatabase(VPstring, [az el], false, 'SOFA');
        H_sum.Data.IR = ifft(fft(H_sum.Data.IR, [], 3) ./ fft(CTF.Data.IR, [], 3), [], 3, 'symmetric');

        % sum the signals
        H_sum.Data.IR = (H_sum.Data.IR + H_sum.Data.IR(1,:,:)) / 2;

        clear pol az el

        %% run the model

        % Sensitivity for localization model (from [1])
        sens = [.21 .26 .44 .44 .46 .52 .55 .58 .63 .70 .71 0.76 0.76 .76 .79 .84 .86 .88 .88 .97 .98 1.02 1.05];
        sens = prctile(sens, 0);

        % model localization performance
        % ('noregular' flag makes sure that nothing is interpolated)
        [pe_qe, pred] = baumgartner2014(H_sum, H_ref, 'QE_PE_EB', 'noregular', 'S', sens, 'lat', lat_ref(ll), 'fs', H_ref.Data.SamplingRate, 'fsstim', H_ref.Data.SamplingRate);

        %% Sort output:
        id = find(round(pred.tang) == round(pol_ref(pp)));
        predp_sorted = [pred.p(:,id:end), pred.p(:,1:id-1)];

        % write to struct:
        baumgartner(:,:,pp,ll) = predp_sorted;
        
        %% plot the results
        if doPlot
            AKf(20,10)
            % normalized localization propability
            imagesc(pol_delta, pred.rang, predp_sorted ./ max(predp_sorted))
            AKpSetColormapAndColorbar('Reds', 'EastOutside', .01, [0 1])

            % plot actual source positions
            hold on
            plot(pol_delta, pol_ref(pp)*ones(size(pol_delta)), 'color', [AKcolors('l') .7], 'LineWidth', 4)
            plot(pol_delta, mod(pol_ref(pp)*ones(size(pol_delta))+pol_delta+90, 360)-90, 'color', [AKcolors('l') .7], 'LineWidth', 4)

            title(sprintf('Reference angle: %i lat, %i pol', lat_ref(ll), pol_ref(pp)))
            ylabel 'Polar angle'
            xlabel '\Delta polar angle'

            set(gca,'YDir','normal')
        end
    end
end

data.p = baumgartner;

end

