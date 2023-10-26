% [r, t, isAudible, tEcho_lat, tMask_lat] = detectReflections(DOA, RIR, fs, timeMax, timeRange, angleRange, thEcho_lat, thEcho_pol, thMask_lat, thMask_pol, eEcho, eMask, t_mix, doPlot, ISM_order, bypassThresold)
% detects potentially audible reflections from a spatial impulse response
% given by the pressure response RIR and the directional component DOA
%
% I N P U T
% DOA          - Direction of arrival vector in carthesion coordinates of
%                size [N x 3] (N: number of samples). X, y, z components
%                are in first, second, and third column, respectively.
% RIR          - Room impulse response of size [N x 1]
%                (N: number of samples)
% fs           - sampling frequency in Hz
% timeMax      - maximum time in ms up to which reflections are identified
%                after the direct sound.
% timeRange    - [1 x 2] time in ms before and after a detected reflection for
%                estimating amplitude and direction
% angleRange   - Determines how many neighboring sound field contributions
%                are considered to form one sound event (direct sound or
%                reflection). This three element vector gives the method in
%                the first element and the parameters for first sound/echos
%                and reflections in the second and third elements
%                [1 w_lat1 w_lat2 w_pol1 w_pol2] : 
%                            uses lateralThreshold.m and polarThreshold.m
%                            to get the lateral and polar range and scales
%                            them with width (w). All contributions within
%                            these ranges are considered for one sound
%                            event. Contributions outside at least one of
%                            the ranges are considered a separate sound
%                            event. A width of 6 will include all
%                            distributions within timeRange.
%                [2 g1 g2] : maximum great circle distance (g) in degree
%                            between main contribution (largest peak) and
%                            any other contribution. A great circle
%                            distance of 180 will include all distributions
%                            within timeRange.
% thEcho_lat   - parameters for defining the lateral echo threshold
%                [slope offset width depth transition_l transition_t1 transition_t2]
%                transition_l smoothness of the transition in lateral
%                direction
%                transition_t1 time in ms where the depth is set to 0 dB
%                transition_t2 time in ms after t1 until the full depth is
%                reached
% thEcho_pol   - parameters for defining the polar echo threshold
%                [slope offset width depth transition_l transition_t1 transition_t2]
% thMask_lat   - parameters for defining the lateral masking threshold
%                [slope offset width depth transition_l transition_t1 transition_t2]
% thMask_pol   - parameters for defining the polar masking threshold
%                [slope offset width depth transition_l transition_t1 transition_t2]
% eEcho        - Determines how the energy of reflections affects the echo
%                threshold
%                [1 scale]: energy is added to the threshold after being
%                           scaled with scale
% eMask        - Determines how the energy of reflections affects the
%                masking threshold (see eEcho for description).
% t_mix        - estimated perceptual mixing time in ms of the room (Is used in
%                the plots.
% doPlot       - false to omitt plots, true or string that specifies the
%                plot title to do the plots
% ISM_order    - for plotting: Order of reflections from the image source
%                model in size [N x 3] (N: Number of ISM reflections) 
%                [t, N_ISM,n] (t: time of reflection in seconds, N: Image
%                source order of reflection, n: sample number of reflection
% bypassThreshold - set true to effectively bypass the masking threshold by
%                   setting the value of the thresholds to 35dB below the
%                   direct sound and the echo threshold to the direct sound
%                
%
% O U T P U T
% r            - a struct containing information about each audible
%                reflection/direct sound contribution:
%       t      - [M x 1] vector holding the onset times of detected direct
%                 sound and reflections in seconds
%       n      - same as t but time in samples
%       a      - [M x 1] vector with the amplitudes of the reflections
%       xzy    - [M x 3] vector with the direction of the reflection in
%                carthesian coordinates
%       ae     - [M x 2] as xyz, but with directions in vertical polar
%                coordinates in degrees [azimuth elevation]
%       lp     - same as ae but with coordinates in interaural polar form
%                [lateralAngle polarAngle]
%       id     - [M x 1] struct with each entry listing the samples that
%                belong to an audible reflection
%       L_echo - gives the level of the reflection relative to the echo
%                threshold in dB. L_echo>1 indicates an audible echo.
%       L_mask - gives the level of the reflection relative to the masking
%                threshold in dB. L_mask>0 indicates an audible reflection.
% t            - list of all processed reflections/contributions similar to
%                r (see above).
% isAudible    - [N x 1] vector that equals 1 of the sample is part of an
%                audible reflection and is 0 otherwise
% tEcho_lat    - [N x 181] matrix holding the temopral evolution of the
%                lateral echo threshold for the lateral angles -90:1:90
% tMask_lat    - [N x 181] matrix holding the temopral evolution of the
%                lateral masking threshold for the lateral angles -90:1:90
% tEcho_pol    - [N x 360] matrix holding the temopral evolution of the
%                polar echo threshold for the polar angles -90:1:269
% tMask_pol    - [N x 360] matrix holding the temopral evolution of the
%                polar masking threshold for the polar angles -90:1:269
%
% written by Tobias Jüterbock, Audio Communication Group TU Berlin, 
% adapted from
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin &
% Microsoft Research, Redmond, USA

%   Copyright 2019 Microsoft Corporation
%   
%   Permission is hereby granted, free of charge, to any person obtaining a 
%   copy of this software and associated documentation files (the "Software"), 
%   to deal in the Software without restriction, including without limitation 
%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%   and/or sell copies of the Software, and to permit persons to whom the 
%   Software is furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in 
%   all copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
%   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
%   DEALINGS IN THE SOFTWARE.

% TRY:
% - When reflections are selected the reference for the angleRange is taken
%   by the single largest peak in the timeRange. This might be changed to
%   the mean directions of all peaks within x dB from the largest peak to
%   be more robust -> variable: levelRange. This would affect the reference
%   DOA and time that would be averaged VBAP style, e.g., rir.^2 * DOA, with
%   sum(rir.^2) == 1
% - The most time consuming part is calculating the time, level, and
%   direction for each sample in the large for loop. A carefully tuned
%   threshold based on the level only of that sample might be used to speed
%   things up.
function [r, t, isAudible, tEcho_lat, tMask_lat, tEcho_pol, tMask_pol] = detectReflections(DOA, RIR, fs, timeMax, timeRange, angleRange, thEcho_lat, thEcho_pol, thMask_lat, thMask_pol, eEcho, eMask, t_mix, doPlot, ISM_order, bypassThresold)

% --------------------------------------------------- input parameter check
if ischar(doPlot)
    plotTitle = doPlot;
    doPlot    = true;
end
if ~exist('plotTitle', 'var')
    plotTitle = '';
end
if ~exist('ISM_order', 'var')
    ISM_order = double.empty(0,3);
end
if ~exist('bypassThresold', 'var')
    bypassThresold = false;
end

thEcho_lat(1:2) = -abs(thEcho_lat(1:2));
thMask_lat(1:2) = -abs(thMask_lat(1:2));
thEcho_pol(1:2) = -abs(thEcho_pol(1:2));
thMask_pol(1:2) = -abs(thMask_pol(1:2));

id              = isnan(sum(DOA, 2));
DOA(id,:)       = repmat([1 0 0], [sum(id) 1]);
RIR(isnan(RIR)) = 0;

timeRange = round(timeRange/1e3*fs); % time range in samples
RIRsq     = RIR.^2;                  % squared IR for calucalting the cumulative energy

% ----------------------------------------------------- truncate to timeMax
N     = min(size(DOA,1), ceil(timeMax/1e3*fs) + timeRange(2) + 1);
DOA   = DOA(1:N,:);
RIR   = RIR(1:N);
ISM_order = ISM_order(ISM_order(:,3) < N,:); %TODO - in large room this truncates one 1st order reflection
clear id

% ----------------------------------------------- allocate space for output

% list of audible contributions
r.t      = nan(N, 1);   % arrival times in seconds
r.n      = nan(N, 1);   % arrival times in samples
r.a      = nan(N, 1);   % amplitudes
r.xyz    = nan(N, 3);   % doa in carthesian coordinates
r.ae     = nan(N, 2);   % doa in vertical polar coordinates [azimuth elevation]
r.lp     = nan(N, 2);   % doa in interaural polar coordinates [later polar]
r.id     = cell(N,1);   % list of samples belonging to one sound event
r.L_echo = nan(N, 1);   % gives reflection the level relative to the echo threshold
r.L_mask = nan(N, 1);   % gives reflection the level relative to the mask threshold
r.TH_mask= nan(N, 1);   % gives the masking threshold value
r.TH_mask_lat = nan(N, 1);   % gives the lateral masking threshold value
r.TH_mask_pol = nan(N, 1);   % gives the polar masking threshold value
% list of all contributions
t.a      = zeros(N, 1); % amplitudes
t.rms    = t.a;         % rms vector
t.xyz    = nan(N, 3);   % doa in carthesian coordinates
t.ae     = nan(N, 2);   % doa in vertical polar coordinates [azimuth elevation]
t.lp     = nan(N, 2);   % doa in interaural polar coordinates [lateral polar]
t.id     = cell(N,1);

% list of reference points
ref.id  = nan(N, 1);
ref.lat = nan(N, 1);
ref.pol = nan(N, 1);
ref.a   = nan(N, 1);
ref.n   = nan(N, 1);

% remaining
isAudible = false(N,1);

% ------------------------------------------------------ detect first sound
disp('detecting first sound')
% get DOA in interaural polar coordinates
[AZ, EL]      = cart2sph(DOA(:,1),DOA(:,2),DOA(:,3));
[LAT, POL]    = sph2hor(AZ./pi.*180, EL./pi.*180);
LP            = [LAT POL];
clear AZ EL LAT POL

% onset
r.n(1) = onsetDetectNR(RIR, fs, fs/4, -30);
r.n(1) = round(r.n(1));

% doa, level, and refined onset
if angleRange(1) == 1
    [r.a(1), r.xyz(1,:), r.ae(1,:), r.lp(1,:), r.n(1), targetReflect] = getDOAandLevel_lp (LP,  RIR, r.n(1), timeRange, angleRange([1 2 4]), isAudible, N); % average [lat pol]
elseif angleRange(1) == 2
    [r.a(1), r.xyz(1,:), r.ae(1,:), r.lp(1,:), r.n(1), targetReflect] = getDOAandLevel_doa(DOA, RIR, r.n(1), timeRange, angleRange([1 2]), isAudible, N); % average [x y z]
else
    error('implementation missing')
end
r.t(1)                   = r.n(1)/fs;
r.L_echo(1)              = abs(thEcho_lat(2));
r.L_mask(1)              = abs(thMask_lat(2));
r.id{1}                  = targetReflect;
isAudible(targetReflect) = true;

% write to target vectors
t.a(r.n(1))     = r.a(1);
t.xyz(r.n(1),:) = r.xyz(1,:);
t.ae(r.n(1),:)  = r.ae(1,:);
t.lp(r.n(1),:)  = r.lp(1,:);
t.id{r.n(1)}    = r.id{1};

% reference id, direction, level, and time for the threshold curve
ref.id(1)  = 1;
ref.lat(1) = r.lp(1,1);
ref.pol(1) = r.lp(1,2);
ref.a(1)   = r.a(1);
ref.n(1)   = r.n(1);

% --------------------------------------- get list of potential reflections
% NOTE: For simplicity the detection of audible reflections is implemented
% in two loops. The first loop (while loop) detects all possible
% reflections and considers large energetic contributions first. The second
% loop (for loop) detects audible reflections from the list generated in
% the first loop.
% Implementing the algorithm in a single loop should not change too much,
% but it has to be made sure that the window used to define a potential
% reflection (using timeRange) is always centered around the largest
% energetic contribution.
target                = false(N, 1);
target(1:r.n(1))      = true;
target(targetReflect) = true;
disp('getting possible reflections');tic
while ~all(target)
    
    id       = find(target == 0);
    [~, id2] = max(abs(RIR(id)));
    n        = id(id2);
    % ---- get level and direction for potential reflection belonging to the current sample
    if angleRange(1) == 1
        [t.a(n), t.xyz(n,:), t.ae(n,:), t.lp(n,:), ~, targetReflect] = getDOAandLevel_lp (LP,  RIR, n, timeRange, angleRange([1 3 5]), target, N); % average [lat pol]
    elseif angleRange(1) == 2
        [t.a(n), t.xyz(n,:), t.ae(n,:), t.lp(n,:), ~, targetReflect] = getDOAandLevel_doa(DOA, RIR, n, timeRange, angleRange([1 3]), target, N); % average [x y z]
    else
        error('implementation missing')
    end
    
    target(targetReflect) = true;
    t.id{n}               = targetReflect;
end
toc
clear target id id2 n

%% ---- blockwise rms --- %
% % check blocksize
% Nrms = 2^10;
% if ~mod(Nrms, 2)
%     Nrms = Nrms+1;
% end
% Nrms2 = floor(Nrms/2);
% 
% % prepare the RIR for blockwise processing
% RIRsq = t.a.^2;
% RIRsq = [repmat(RIRsq(1,:), Nrms2, 1); RIRsq; repmat(RIRsq(end,:), Nrms2+1, 1)];
% 
% % initialize the sum of squares
% sq       = sum(RIRsq(1:Nrms,:));
% 
% % blocwise estimate
% for nn = 1:size(RIRsq,1)-Nrms
%     
%     t.rms(nn,:) = sq;
%     sq = sq - RIRsq(nn,:) + RIRsq(nn+Nrms,:);
%     
% end
% 
% t.rms = sqrt(t.rms / Nrms);
% 
% clear Nrms Nrms2 RIRsq sq nn

%% ---------------------------------------------- detect audible reflections

% counters
nReflections = 1;
nEchos       = 1;

% convert timeMax to samples
nMax = min(r.n(1)+round(timeMax/1e3*fs)+timeRange(2), N);

% allocate space
tEcho_lat = nan(N, 181);
tEcho_pol = nan(N, 360);
tMask_lat = tEcho_lat;
tMask_pol = tEcho_pol;
tMask     = nan(N,1);

cumE  = nan(N, 1);

% get the first entry of the masking and echo threshold
tEcho_lat(ref.n(nEchos),:) = ref.a(nEchos) * 10^( thEcho_lat(2) / 20 );
tEcho_pol(ref.n(nEchos),:) = ref.a(nEchos) * 10^( thEcho_pol(2) / 20 );
tMask_lat(ref.n(nEchos),:) = ref.a(nEchos) * 10^( thMask_lat(2) / 20 );
tMask_pol(ref.n(nEchos),:) = ref.a(nEchos) * 10^( thMask_pol(2) / 20 );

disp('detecting audible reflections')
for nn = ref.n(nEchos)+1:nMax
    
    % get cumulative energy from the latest previous reference up to the current potential reflection
    % NOTE: with the current implementation the current reference point
    % given by nEcho is always the one that is used. If the algorithm is
    % implemented in a single loop, this might not always be the case.
    id                            = find(ref.n < nn, 1, 'last');
    idE                           = false(size(RIRsq));
    idE(min(r.id{ref.id(id)}):nn) = 1; % Include bins between reference (first sound or echo) and current sample
    idE(r.id{ref.id(id)})         = 0; % exclude energy from reference
    
    cumE(nn) = sqrt(sum(RIRsq(idE)));
    
    if bypassThresold == false
    % thresholds for all angles 
    % TODO: tEcho_polar und tMask_polar auch weiter verarbeiten
    tEcho_lat(nn,:) = getTH(ref.n(id), nn, ref.lat(id),ref.pol(id), -90:90,  ref.pol(id), ref.a(id), thEcho_lat, thEcho_pol, cumE(nn), eEcho, fs, 'lateral'); 
    tEcho_pol(nn,:) = getTH(ref.n(id), nn, ref.lat(id),ref.pol(id), ref.lat(id), -90:269, ref.a(id), thEcho_lat, thEcho_pol, cumE(nn), eEcho, fs, 'polar'  ); 
    tMask_lat(nn,:) = getTH(ref.n(id), nn, ref.lat(id),ref.pol(id), -90:90,  ref.pol(id), ref.a(id), thMask_lat, thMask_pol, cumE(nn), eMask, fs, 'lateral');
    tMask_pol(nn,:) = getTH(ref.n(id), nn, ref.lat(id),ref.pol(id), ref.lat(id), -90:269, ref.a(id), thMask_lat, thMask_pol, cumE(nn), eMask, fs, 'polar'  );
    end
    % make sure the masking threshold does not exceed the echo threshold
    thRatio_lat = min(tEcho_lat(nn,:)) / min(tMask_lat(nn,:));
    thRatio_pol = min(tEcho_pol(nn,:)) / min(tMask_pol(nn,:));
    if thRatio_lat < 1
        tMask_lat(nn,:) = tMask_lat(nn,:) * thRatio_lat;
    end
    if thRatio_pol < 1
        tMask_pol(nn,:) = tMask_pol(nn,:) * thRatio_pol;
    end
    
    % check if the current sample is already part of r.a reflection && if
    % it is an onset of a potential reflection
    if ~isAudible(nn) && ~isempty(t.id{nn})
        
        % get thresholds
        if bypassThresold == false
        echoTH     = getTH(ref.n(id), nn, ref.lat(id), ref.pol(id), t.lp(nn,1),t.lp(nn,2), ref.a(id), thEcho_lat, thEcho_pol, cumE(nn), eEcho, fs, 'both');
%         echoTH_lat = getTH(ref.n(id), nn, ref.lat(id), ref.pol(id), t.lp(nn,1),t.lp(nn,2), ref.a(id), thEcho_lat, thEcho_pol, cumE(nn), eEcho, fs, 'lateral');
%         echoTH_pol = getTH(ref.n(id), nn, ref.lat(id), ref.pol(id), t.lp(nn,1),t.lp(nn,2), ref.a(id), thEcho_lat, thEcho_pol, cumE(nn), eEcho, fs, 'polar');
        maskTH     = getTH(ref.n(id), nn, ref.lat(id), ref.pol(id), t.lp(nn,1),t.lp(nn,2), ref.a(id), thMask_lat, thMask_pol, cumE(nn), eMask, fs, 'both');
        maskTH_lat = getTH(ref.n(id), nn, ref.lat(id), ref.pol(id), t.lp(nn,1),t.lp(nn,2), ref.a(id), thMask_lat, thMask_pol, cumE(nn), eMask, fs, 'lateral');
        maskTH_pol = getTH(ref.n(id), nn, ref.lat(id), ref.pol(id), t.lp(nn,1),t.lp(nn,2), ref.a(id), thMask_lat, thMask_pol, cumE(nn), eMask, fs, 'polar');
        tMask(nn) = maskTH;
        if abs(maskTH_lat - maskTH) > 1e-13 && abs(maskTH_pol - maskTH) > 1e-13
            fprintf('Lateral threshold: %f \n',db(maskTH_lat))
            fprintf('polar threshold: %f \n',db(maskTH_pol))
            fprintf('Threshold: %f \n',db(maskTH))
            fprintf(2,'Threshold not equal to lat or pol! \n\n\n')
            pause(2)
        end
        elseif bypassThresold == true
            % set masking threshold to 35dB below direct sound
            maskTH =  10.^((db(r.a(1))-35) / 20);
            maskTH_lat = maskTH;
            maskTH_pol = maskTH;
            echoTH = 10.^(db(r.a(1)) / 20); 
        end

        % make sure the masking threshold does not exceed the echo threshold
        if thRatio_lat < 1
            maskTH = maskTH * thRatio_lat;
        end
        
        % test thresholds
        if t.a(nn) > echoTH
            % set write flags
            writeReflection = true;
            writeEcho       = true;
        else
            writeEcho = false;
            
            if t.a(nn) > maskTH
                % set write flag
                writeReflection = true;
            else
                writeReflection = false;
            end
        end
        
        % write reflection data
        if writeReflection
            nReflections = nReflections + 1;
            
            r.a(nReflections)      = t.a(nn);
            r.n(nReflections)      = nn;
            r.t(nReflections)      = nn/fs;
            r.xyz(nReflections,:)  = t.xyz(nn,:);
            r.ae(nReflections,:)   = t.ae(nn,:);
            r.lp(nReflections,:)   = t.lp(nn,:);
            r.id{nReflections}     = t.id{nn};
            r.L_echo(nReflections) = db(t.a(nn)/echoTH);
            r.L_mask(nReflections) = db(t.a(nn)/maskTH);
            r.TH_mask(nReflections)= db(maskTH); 
            r.TH_mask_lat(nReflections)= db(maskTH_lat); 
            r.TH_mask_pol(nReflections)= db(maskTH_pol); 
            
            isAudible(t.id{nn})  = 1;
        end
        % write echo data
        if writeEcho
            nEchos = nEchos + 1;
            
            % set new reference
            ref.id(nEchos)  = nReflections;
            ref.lat(nEchos) = r.lp(nReflections);
            ref.a(nEchos)   = r.a(nReflections);
            ref.n(nEchos)   = r.n(nReflections);
        end   
    end 
end

% ----------------------------------------------------------- format output
r.t       = r.t  (1:nReflections);
r.n       = r.n  (1:nReflections);
r.a       = r.a  (1:nReflections);
r.xyz     = r.xyz(1:nReflections,:);
r.ae      = r.ae (1:nReflections,:);
r.lp      = r.lp (1:nReflections,:);
r.id      = r.id(1:nReflections);
r.L_echo  = r.L_echo(1:nReflections);
r.L_mask  = r.L_mask(1:nReflections);
r.TH_mask = r.TH_mask(1:nReflections);
r.TH_mask_lat = r.TH_mask_lat(1:nReflections);
r.TH_mask_pol = r.TH_mask_pol(1:nReflections); % TODO: remove this when done with debugging

% -------------------------------------------------------------------- plot
if doPlot
    plotDetectedReflections(DOA, RIR, isAudible, fs, r, t, t_mix, N, tEcho_lat, tEcho_pol, tMask_lat, tMask_pol, timeMax, ISM_order,tMask)
    plot_th(r, t, tEcho_lat, tEcho_pol, tMask_lat, tMask_pol, fs, N, t_mix, plotTitle);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, xyz, ae, lp, timeSample, reflections] = getDOAandLevel_lp(LP, RIR, timeSample, timeRange, angleRange, isAudible, N)
    
% --------------------------- get samples that contribute to the reflection
% ---- select according to timeRange
id  = max(timeSample-timeRange(1),1):min(timeSample+timeRange(2),N);
id2 = find(isAudible(id) == 0);
lp  = LP(id(id2),:);
rir = RIR(id(id2));

% ---- select according to angleRange
% get reference direction
[az, el] = hor2sph(lp(:,1), lp(:,2));
[~, id3] = max(abs(rir));

if     angleRange(1) == 1
    % select based on lateral and polar angle
    [~, latRef] = lateralThreshold([], lp(id3,1), angleRange(2));
%     tic
    [~, polRef] = polarThreshold([], lp(id3,1),lp(id3,2), angleRange(3));
%     toc
    if     polRef(1) > polRef(2)
        id4     = lp(:,1)<=latRef(1) & lp(:,1)>=latRef(2) & ...
                  lp(:,2)<=polRef(1) & lp(:,2)>=polRef(2);
    elseif polRef(1) < polRef(2)
       id4     = (lp(:,1)<=latRef(1) & lp(:,1)>=latRef(2)) &      ...
                 (lp(:,2)>=polRef(2) | lp(:,2)<=polRef(1));
    elseif polRef(1) == polRef(2)
        error(['polRef upper and lower limit are equal. This should only',...
              'happen if angleRange is set to 0.'])
    end

elseif angleRange(1) == 2
    % select based on great circle distance (gcd)
    gcd = acos( sin(el)*sin(el(id3)) +  cos(el)*cos(el(id3)).*cos(az-az(id3)) ) / pi * 180;
    id4 = find(gcd <= angleRange(2));
else
    error('No such method')
end
lp  = lp(id4,:);
rir = rir(id4);

% ------------- return the position of the maximum value as the final onset
timeSample = id( id2( id3 ) );

% -------------------------------- mark audible contributions for returning
reflections = id( id2( id4 ) );

% ------------------------------------------------------- get the amplitude
a = sqrt(sum(rir.^2));

% --------------------- get the doa of the reflection by weighted averaging

% get doa VBAP style and using a weighted circular mean for the polar angle
lat = rir.^2' * lp(:,1) / sum(rir.^2);
pol = angle( rir.^2' * exp( 1j * lp(:,2)/180*pi ) ) / pi*180;
pol = mod(pol+90,360)-90;

lp = [lat pol];

% -------------------------------------------- convert to other coordinates
[az, el]  = hor2sph(lat, pol);
ae        = [az el];

[x, y, z] = sph2cart(az/180*pi, el/180*pi, 1);
xyz       = [x y z];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, xyz, ae, lp, timeSample, reflections] = getDOAandLevel_doa(DOA, RIR, timeSample, timeRange, angleRange, isAudible, N)
    
% --------------------------- get samples that contribute to the reflection
% ---- select according to timeRange
id  = max(timeSample-timeRange(1),0):min(timeSample+timeRange(2),N);
id2 = find(isAudible(id) == 0);
doa = DOA(id(id2),:);
rir = RIR(id(id2));

% ---- select according to angleRange
% get reference direction
[az, el] = cart2sph(doa(:,1), doa(:,2), doa(:,3));
[~, id3] = max(abs(rir));

if     angleRange(1) == 1
    % select based on lateral angle
    [lat, pol]  = sph2hor(az/pi*180, el/pi*180);
    [~, latRef] = lateralThreshold([], lat(id3), angleRange(2));
    [~, polRef] = polarThreshold([], lat(id3), pol(id3), angleRange(4));
    if     polRef(1) > polRef(2)
        id4     = (lat<=latRef(1) & lat>=latRef(2)) & ...
                  (pol<=polRef(1) & pol>=polRef(2));
    elseif polRef(1) < polRef(2)
        id4     = (lat<=latRef(1) & lat>=latRef(2)) & ...
                  (pol>=polRef(1) | pol<=polRef(2));
    elseif polRef(1) == polRef(2)
        error(['polRef upper and lower limit are equal. This should only',...
              'happen if angleRange is set to 0.'])
    end
elseif angleRange(1) == 2
    % select based on great circle distance (gcd)
    gcd = acos( sin(el)*sin(el(id3)) +  cos(el)*cos(el(id3)).*cos(az-az(id3)) ) / pi * 180;
    id4 = find(gcd <= angleRange(2));
else
    error('No such method')
end
doa      = doa(id4,:);
rir      = rir(id4);

% ------------- return the position of the maximum value as the final onset
timeSample = id( id2( id3 ) );

% -------------------------------- mark audible contributions for returning
reflections = id( id2( id4 ) );

% ------------------------------------------------------- get the amplitude
a = sqrt(sum(rir.^2));

% --------------------- get the doa of the reflection by weighted averaging

% get doa VBAP style
xyz = rir.^2' * doa;
xyz = xyz / sqrt(sum(xyz.^2));

% ---------------------------------------- convert to spherical coordinates
[az, el]  = cart2sph(xyz(:,1), xyz(:,2), xyz(:,3));
ae        = [az el] / pi * 180;
ae(:,1)   = mod(ae(:,1), 360);

[lat, pol] = sph2hor(ae(:,1), ae(:,2));
lp         = [lat, pol];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function th = getTH(refN, targetN, refLat, refPol, latTarget, polTarget, refA, TH_lat, TH_pol, cumEnergy, energyMethod, fs, thresholdMethod)
%NOTE: latTarget and polTarget have expected to be pairs of target angles
%(they need to have the same number of elements)

% time in ms since first sound / echoed sound
delta_t = (targetN - refN) / fs * 1e3;

% ---- get the current echo threshold
% depth
if delta_t <= TH_lat(6)
    depth_lat = 0;
elseif delta_t <= TH_lat(6)+TH_lat(7)
    depth_lat = TH_lat(4) * (delta_t-TH_lat(6)) / TH_lat(7);
else
    depth_lat = TH_lat(4);
end
if delta_t <= TH_pol(6)
    depth_pol = 0;
elseif delta_t <= TH_pol(6)+TH_pol(7)
    depth_pol = TH_pol(4) * (delta_t-TH_pol(6)) / TH_pol(7);
else
    depth_pol = TH_pol(4);
end

% threshold at current time
if energyMethod(1) == 1
    if strcmp(thresholdMethod, 'lateral')
        % lateral threshold in dB                       
        th  = lateralThreshold(latTarget, refLat, TH_lat(3), depth_lat, TH_lat(5));
        % linear threshold
        th  = 10.^(th / 20);
        % amplitude: reference * time dependency * offset
        thA = refA * 10^( delta_t*TH_lat(1) / 20 ) * 10^( TH_lat(2) / 20 );
        % add cumulative energy
        thA = (thA + cumEnergy*energyMethod(2));
        % apply amplitude
        th  = th * thA;
        
    elseif strcmp(thresholdMethod, 'polar')
        % polar threshold in dB                                
        th  = polarThreshold(polTarget, refLat, refPol, TH_pol(3), depth_pol, TH_pol(5)); % TODO was ist TH(3) und (5)
        % linear threshold
        th  = 10.^(th / 20);
        % amplitude: reference * time dependency * offset
        thA = refA * 10^( delta_t*TH_pol(1) / 20 ) * 10^( TH_pol(2) / 20 );
        % add cumulative energy
        thA = (thA + cumEnergy*energyMethod(2));
        % apply amplitude
        th  = th * thA;
    
    elseif strcmp(thresholdMethod, 'both')
        % lateral threshold in dB                       
        th  = lateralThreshold(latTarget, refLat, TH_lat(3), depth_lat, TH_lat(5));
        % linear threshold
        th  = 10.^(th / 20);
        % amplitude: reference * time dependency * offset
        thA = refA * 10^( delta_t*TH_lat(1) / 20 ) * 10^( TH_lat(2) / 20 );
        % add cumulative energy
        thA = (thA + cumEnergy*energyMethod(2));
        % apply amplitude
        th_lat  = th * thA;

        % polar threshold in dB                                
        th  = polarThreshold(polTarget, refLat, refPol, TH_pol(3), depth_pol, TH_pol(5)); % TODO was ist TH(3) und (5)
        % linear threshold
        th  = 10.^(th / 20);
        % amplitude: reference * time dependency * offset
        thA = refA * 10^( delta_t*TH_pol(1) / 20 ) * 10^( TH_pol(2) / 20 );
        % add cumulative energy
        thA = (thA + cumEnergy*energyMethod(2));
        % apply amplitude
        th_pol  = th * thA;
        % TODO: FIXME: target jeweils für lat und pol machen
        th = min(th_lat,th_pol);
    else
        error('select thresholdMethod ''polar'', ''lateral'' or ''both'' ')
    end
    
else
    error('not defined')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotDetectedReflections(DOA, RIR, isAudible, fs, r, t, t_mix, N, tEcho_lat, tEcho_pol, tMask_lat, tMask_pol, timeMax, ISM_order,tMask)

%% audibility and maximum
idAud = logical(isAudible);
[~, mID] = max(abs(RIR));

% maximum plot length
nMax = min(round(r.t(1)*fs+.1*fs), numel(RIR));

% figure
ph = newFig(12,15);
ph.Name = 'SRIR';

subplot('Position', [0.11,0.7167,0.85,0.2733])
    % Thresholds
    idEcho = ~isnan(tEcho_lat(:,1));
    Echo_lat   = db( [min(tEcho_lat(idEcho,:), [], 2) max(tEcho_lat(idEcho,:), [], 2)] );
    Echo_pol   = db( [min(tEcho_pol(idEcho,:), [], 2) max(tEcho_pol(idEcho,:), [], 2)] );
    Mask_lat   = db( [min(tMask_lat(idEcho,:), [], 2) max(tMask_lat(idEcho,:), [], 2)] );
    Mask_pol   = db( [min(tMask_pol(idEcho,:), [], 2) max(tMask_pol(idEcho,:), [], 2)] );
    % time vector
    time = (0:N-1)'/fs;
    time = time(idEcho);
    % max. for normalization
    aMax = db(max(abs(t.a)));
    % plot threshold ranges
    h1 = patch([time; flipud(time)]*1e3, [Echo_lat(:,1); flipud(Echo_lat(:,2))]-aMax, [0.4940    0.1840    0.5560], 'EdgeColor', [0.4940    0.1840    0.5560], 'FaceAlpha', .5, 'EdgeAlpha', .5, 'LineWidth', 1);
    h2 = patch([time; flipud(time)]*1e3, [Mask_lat(:,1); flipud(Mask_lat(:,2))]-aMax, [.266 .674 .188], 'EdgeColor', [.266 .674 .188], 'FaceAlpha', .3, 'EdgeAlpha', .5, 'LineWidth', 1);
    h3 = patch([time; flipud(time)]*1e3, [Mask_pol(:,1); flipud(Mask_pol(:,2))]-aMax, [.850,.325,.098], 'EdgeColor', [.850,.325,.098], 'FaceAlpha', .3, 'EdgeAlpha', .5, 'LineWidth', 1);
    hold on
    % plot inaudible
    rir     = RIR/max(abs(RIR));
    rir(idAud) = 0;
    rir(mID) = RIR(mID)/max(abs(RIR))*.95; % copy the maximum to keep the scale and scale to make sure that it is not plotted above the actual maximum
    plotSRIR(rir, DOA, 'rir', fs, r.t(1)+.1, 30, false, 'scatter', {'.' 1 10}, 'B');
    % plot audible
    rir      = RIR/max(abs(RIR));
    rir(~idAud) = 0;
    plotET(rir(1:nMax), fs, false, [.94 .23 .17], '-', false, [-30 5]);
    % instead of: AKp(rir(1:round(r.t(1)*fs+.1*fs)), 'et2d', 'dr', [-30 5], 'fs', fs, 'c', [.94 .23 .17])
    % mark level and echo
    hold on
    plot(r.t*1e3, db(r.a/max(abs(r.a))), 'ko', 'LineWidth', 1.2)
    for nn = 1:numel(r.L_echo)
        if r.L_echo(nn)>0
            plot(r.t(nn)*1e3, db(r.a(nn)/max(abs(r.a))), 'ko', 'LineWidth', 1.2, 'Color', [0.85 0.225 0.098])
        end
    end
    % plot mixing time
    plot([t_mix t_mix], [-30 5], '--k', 'Linewidth', 1.2, 'color', [.6 .6 .6],'HandleVisibility','off')
    
    xlabel ''; title ''
    set(gca, 'XTickLabel', '')
    xlim([0 timeMax])
    legend([h1 h2 h3], 'Echo threshold', 'lat. mask. threshold','pol. mask. threshold')
    
subplot('Position', [0.11,0.3933,0.85,0.2733])
    % plot inaudible
    rir      = RIR/max(abs(RIR));
    rir(idAud)  = 0;
    rir(mID) = RIR(mID)/max(abs(RIR))*.95; % copy the maximum to keep the scale and scale to make sure that it is not plotted above the actual maximum
    h1 = plotSRIR(rir, DOA, 'lat', fs, r.t(1)+.1, 30, false, 'scatter', {'.' 1 10}, 'B');
    % plot audible
    rir      = RIR/max(abs(RIR));
    rir(~idAud) = 0;
    h2 = plotSRIR(rir, DOA, 'lat', fs, r.t(1)+.1, 30, false, 'scatter', {'.' 1 10}, 'R', false);
    view([0 90])
    % mark reflection angle
    hold on
    plot3(r.t*1e3, r.lp(:,1), repmat(max(get(gca, 'ZLim')), numel(r.t), 1), 'ko', 'LineWidth', 1.2)
    % plot mixing time
    plot([t_mix t_mix], [-90 90], '--k', 'Linewidth', 1.2, 'color', [.6 .6 .6])
    % plot 1st order reflections
    if ~isempty(ISM_order)
        rir      = zeros(size(rir));
        idx = ISM_order(:,2) == 1;
        idx2     = ISM_order(idx,3);
        rir(idx2)= 1;
        h3 = plotSRIR(rir, DOA, 'lat', fs, r.t(1)+.1, 30, false, 'scatter', {'o' 1 10}, 'G', false);
    end
    xlabel ''
    set(gca, 'XTickLabel', '')    
    xlim([0 timeMax])
    if ~isempty(ISM_order)
        legend([h2 h1 h3],'audible','inaudible','1st order')
    else
        legend([h2 h1],'audible','inaudible')
    end

subplot('Position', [0.11,0.07,0.85,0.2733])
    % plot inaudible
    rir      = RIR/max(abs(RIR));
    rir(idAud)  = 0;
    rir(mID) = RIR(mID)/max(abs(RIR))*.95; % copy the maximum to keep the scale and scale to make sure that it is not plotted above the actual maximum
    plotSRIR(rir, DOA, 'pol', fs, r.t(1)+.1, 30, false, 'scatter', {'.' 1 10}, 'B');
    % plot audible
    rir      = RIR/max(abs(RIR));
    rir(~idAud) = 0;
    plotSRIR(rir, DOA, 'pol', fs, r.t(1)+.1, 30, false, 'scatter', {'.' 1 10}, 'R', false);
    view([0 90])
    % mark reflection angle
    hold on
    plot3(r.t*1e3, r.lp(:,2), repmat(max(get(gca, 'ZLim')), numel(r.t), 1), 'ko', 'LineWidth', 1.2)
    % plot mixing time
    plot([t_mix t_mix], [-90 270], '--k', 'Linewidth', 1.2, 'color', [.6 .6 .6])
    % plot 1st order reflections
    if ~isempty(ISM_order)
        rir      = zeros(size(rir));
        idx = ISM_order(:,2) == 1;
        idx2     = ISM_order(idx,3);
        rir(idx2)= 1;
        plotSRIR(rir, DOA, 'pol', fs, r.t(1)+.1, 30, false, 'scatter', {'o' 1 10}, 'G', false);
    end
    xlim([0 timeMax])
    xTicks = xticks;
    xticks(xTicks(1:end-1))
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ph = plot_th(r, t, tEcho_lat, tEcho_pol, tMask_lat, tMask_pol, fs, N, t_mix, plotTitle)
%%
% echo/masking threshold
id   = ~isnan(tEcho_lat(:,1));
Echo_lat = db( [min(tEcho_lat(id,:), [], 2) max(tEcho_lat(id,:), [], 2)] );
Echo_pol = db( [min(tEcho_pol(id,:), [], 2) max(tEcho_pol(id,:), [], 2)] );
Mask_lat = db( [min(tMask_lat(id,:), [], 2) max(tMask_lat(id,:), [], 2)] );
Mask_pol = db( [min(tMask_pol(id,:), [], 2) max(tMask_pol(id,:), [], 2)] );

% time vector
time = (0:N-1)'/fs;
time = time(id);

% max. for normalization
aMax = db(max(abs(t.a)));

ph = newFig(12,12);

% Plot threshold, reflections, and mark detected reflections
subplot(2,1,1)
    plotET(t.a, fs, true, 'k', '-', false, [-60 0])
    % instead of: AKp(t.a, 'et2d', 'fs', fs, 'dr', 60, 'norm_d', 0)
    hold on
    
    % plot threshold Minima and Maxima
    h1 = patch([time; flipud(time)]*1e3, [Echo_lat(:,1); flipud(Echo_lat(:,2))]-aMax, [0.4940    0.1840    0.5560], 'EdgeColor', [0.4940    0.1840    0.5560], 'FaceAlpha', .5, 'EdgeAlpha', .5, 'LineWidth', 1);
%     h2 = patch([time; flipud(time)]*1e3, [Echo_pol(:,1); flipud(Echo_pol(:,2))]-aMax, [0.85,0.325,0.098], 'EdgeColor', [0.85,0.225,.098], 'FaceAlpha', .5, 'EdgeAlpha', .5, 'LineWidth', 1);
    h3 = patch([time; flipud(time)]*1e3, [Mask_lat(:,1); flipud(Mask_lat(:,2))]-aMax, [0    0.447 0.741], 'EdgeColor', [0    0.447 0.741], 'FaceAlpha', .5, 'EdgeAlpha', .5, 'LineWidth', 1);
    h4 = patch([time; flipud(time)]*1e3, [Mask_pol(:,1); flipud(Mask_pol(:,2))]-aMax, [.850,0.325,0.098], 'EdgeColor', [0.85 0.326 0.098], 'FaceAlpha', .5, 'EdgeAlpha', .5, 'LineWidth', 1);
    
    % plot actual Masking thresholds:
    plot(r.t(2:end)*1e3, r.TH_mask(2:end)-aMax,'-','Color',[.1,.6,.2,0.3],'LineWidth',1.2,'MarkerSize',7)
    h5 = plot(r.t(2:end)*1e3, r.TH_mask(2:end)-aMax,'.','Color',[.1,.7,.2,1],'LineWidth',1.2,'MarkerSize',7);
    % plot detected reflections
    plot(r.t*1e3, db(r.a)-aMax, '.', 'color', [0.85 0.225 0.098], 'LineWidth', 1.2, 'MarkerSize', 5)
    
    % plot mixing time
    Lims = [get(gca, 'XLim') get(gca, 'YLim')];
    plot([t_mix t_mix], Lims(3:4), '--k', 'Linewidth', 1.2, 'color', [.6 .6 .6])

    title(['Reflection Thresholds (' plotTitle ')'])

% plot level od detected Reflections relative to masking threshold
subplot(2,1,2)
    % plot reflections
    plot(r.t(2:end)*1e3, r.L_mask(2:end), '.-k', 'LineWidth', 1.2, 'Color', [0.85 0.225 0.098 .5], 'MarkerEdgeColor', [0.85 0.225 0.098], 'MarkerSize', 15)
    xlim(Lims(1:2))
    
    % plot mixing time
    hold on
    Lims = [get(gca, 'XLim') get(gca, 'YLim')];
    plot([t_mix t_mix], Lims(3:4), '--k', 'Linewidth', 1.2, 'color', [.6 .6 .6])
    
    xlabel 'Time in ms'
    ylabel 'Relative amplitude in dB'
    title('reflection levels, relative to masking threshold')
    
    box on
    grid on
    hh = legend([h1,h3,h4,h5],'Echo Threshold','lat. mask. threshold','pol. mask. threshold','abs. masking threshold','location','NorthOutside');
    hh.Location = 'SouthEast';
ph.Name = 'Threshold';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat,pol]=sph2hor(az,el)

[x,y,z] = sph2cart(deg2rad(az),deg2rad(el),ones(size(az)));

% remove noise below eps
x(abs(x)<eps)=0;
y(abs(y)<eps)=0;
z(abs(z)<eps)=0;

% interpret horizontal polar format as rotated spherical coordinates with
% negative azimuth direction
[pol,nlat] = cart2sph(x,z,-y);
pol        = rad2deg(pol);
lat        = rad2deg(-nlat);

% adjust polar angle range
pol = mod(pol+90,360)-90;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [az,el]=hor2sph(lat,pol)

% interpret horizontal polar format as rotated spherical coordinates with
% negative azimuth direction
[x,nz,y] = sph2cart(-deg2rad(pol),deg2rad(lat),ones(size(lat)));

[az,el] = cart2sph(x,y,-nz);

az = rad2deg(az);
el = rad2deg(el);

% adjust azimuth range 
az = mod(az, 360);
end