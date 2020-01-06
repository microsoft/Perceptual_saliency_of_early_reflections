% [BRIR, ER, LR, LRmismatch, LRdecay, LRdecayOrg] = composeParametricBRIR(RIR, DOA, fs, r, isAudible, headRot, N, tailMethod, freqMethod, spatMethod, doPlot)
% returns a prametric BRIR composed from first sound, early, refelction and
% late reverberation. The first sound and early refelctions are composed
% from reflections in input 'r' using the FABIAN HRTF database [1]. The
% late reverberation is computed based on the blockwise residual RMS energy
% of the RIR, i.e., not considering the energy contained in the first sound
% and early refelctions. The resiudal RMS energy is used to compute a decay
% curve that is applied to gaussian noise. The output sampling rate of the
% parametric BRIR is 44.1 kHz and corresponds to the HRTFs used.
%
% I N P U T:
% DOA        - Direction of arrival vector in carthesion coordinates of
%              size [N x 3] (N: number of samples). X, y, and z components
%              are in first, second, and third column, respectively.
% RIR        - Room impulse response of size [N x 1] (N: number of samples)
% fs         - sampling frequency in Hz.
% r          - List of refelctions in the format according to
%              detectreflections.m and reduceReflections.m
% isAudible  - logical vector of size [N x 1] where true indicates that a
%              sample is part of the first sound or an early refelction.
%              Samples flagged with true will be excluded from calculating
%              the energy of the late reverberation.
% headRot    - two element vector to apply head rotation given by
%              [azimuth rotatin in degree x elevation roation in degree]
% N          - processing block length for calculating the residual RMS
%              energy in samles. (Hop size = 1 sample).
% tailMethod - A cell array given the tail method and it's parameterization
%              by {'method' param1 param2 ... paramN}. The following list
%              gives the methods and parameters by number:
%              'exact'  Directly applies the residual RMS curve to the
%                       gaussian noise. No additional parameters. Use
%                       {'exact'} in function call.
%              'ramp_a' Estimates the decay from the logarithmic residual
%                       energy by fitting a plolynomial of order param1.
%                       The estimation start after the the last
%                       reflection specified in 'r'. For example call
%                       {'ramp' 1} to fit a first order polynom. Warning:
%                       residual energ before the last included reflection
%                       is discarded. A second parameter can be passed to
%                       detrmine the length of the late reverberation in
%                       seconds. E.g., pass {'ramp' 1 1} to generate 1 s of
%                       late reverberation. If param2 is omitted, the late
%                       reverberation matches the duration of the RIR
%                       passed to this function.
%              'ramp_b' As 'ramp_a' but applies a generic  exponential fade
%                       in between the first sound and last included
%                       refelction to account for residual energy.
%              'ramp_c' as 'ramp_a' but applies an exponential fade in
%                       between the first sound and last refelction that
%                       matches the residual energy.
%              'ramp_d' As 'ramp_a' but fits a second polynomial to the log
%                       residual energy starting at the first sound and
%                       ending at the last included refelction to
%                       acount for residual energy.
% freqMethod - Specify frequency dependent processing (to be implemented in
%              future versions - pass false for now)
% spatMethod - Specify spatial resolution of the late reverberation (to be
%              implemented in future versions - pass false for now)
% doPlot     - false to ommit plots, true, or plot name as string to show
%              plots.
%
%
% O U T P U T:
% BRIR       - The parametric binaural room impulse response
% ER         - The early refelctions of the parametric BRIR
% LR         - The late reverberation of the parametric BRIR
% LRmismatch - Energetic missmatch between the parametric and actual late
%              reverberation in dB. Mismatch comes from approximating the
%              actual LR with a simple polygon.
% LRdecay    - The decay curve applied to the late reverberation.
% LRdecayOrg - The decay curve resampled to the input sampling rate (using
%              simple linear interpolation.
%
%
% [1] Fabian Brinkmann et al. (2017): The Fabian head-reated transfer
%     function database. http://dx.doi.org/10.14279/depositonce-5718.3
%
%
% WARNINGS:
% - THIS FUNCTION ASSUMES THAT THE RIR DOES NOT CONTAIN NOISE AT THE END.
%   MAKE SURE THAT THIS IS THE CASE TO AVOID ERROR IN ESTIMATING THE SLOPE
%   OF THE LATE REVERBERATION IF USING METHODS 'ramp' OR 'ramps'.
% - THE RANDOM GENERATOR IS NOT SEEDED WITH THIS FNCTION TO ASSURE THAT THE
%   LATE REVERBERATION SOUNDS IDENTICAL ACROSS FUNCTION CALLS RUN
%   "rng(seed)".
% - THIS FUNCTION USES THE FABIAN HRTF DATABASE AND HAS HARD CODED
%   VARIABLES THAT DESCRIBE PROBERTIES OF THAT DATASET. IF YOU WANT TO USE
%   DIFFRERENT HRTFS CHANGE ALL VARIABLES STARTING WITH 'HRTF_' AND REPLACE
%   THE LINE CALLING 'AKhrirInterpolation.m'.
%
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
function [BRIR, ER, LR, LRmismatch, LRdecay, LRdecayOrg] = composeParametricBRIR(RIR, DOA, fs, r, isAudible, headRot, N, tailMethod, freqMethod, spatMethod, doPlot, show_plot)
%% ------------------------------------------------------------ check input
if ischar(doPlot)
    plotTitle = [' - ' doPlot];
    doPlot    = true;
end
if ~exist('plotTitle', 'var')
    plotTitle = '';
end

% separate by early and late reverberation
isAudible(end+1:numel(RIR)) = false;

RIR_LR             = RIR;
RIR_LR(isAudible)  = 0;
RIR_ER             = RIR;
RIR_ER(~isAudible) = 0;

%% ------------------------------------------------ get early reverberation
% parameters for binaural rendering
HRTF_delay   = 28;       % dely in frontal HRTF in samples
HRTF_fs      = 44100;    % sampling rate of HRTFs

% append zeros to match the HRTF delay
ER = zeros(HRTF_delay+1+round(numel(RIR_LR)*HRTF_fs/fs), 2);

% rotate reflections
[r.ae(:,1), r.ae(:,2)] = headRotation(r.ae(:,1), r.ae(:,2), headRot(1), headRot(2));

% get HRTFs - faster outside the loop for a few reflections
[h_l, h_r] = hrirInterpolation(r.ae(:,1), r.ae(:,2), 0, 'measured_sh_dir');
n          = round(r.t*HRTF_fs)+1 - HRTF_delay;

% HRTF length in samples
HRTF_N     = size(h_l, 1);

% write HRTFs to early refelctions
for mm = 1:numel(r.t)
    ER(n(mm):n(mm)+HRTF_N-1,:) = ER(n(mm):n(mm)+HRTF_N-1,:) + [h_l(:,mm) h_r(:,mm)]*r.a(mm);
end

% remove HRTF_delay
ER = ER(HRTF_delay+2:end,:);

%% ------------------------------------------------- get late reverberation

% ---- get the blockwise RMS in the input sample rate fs ---- %
% blockwise rms is designed for odd block lengths
if ~mod(N, 2)
    N = N+1;
end
N2 = floor(N/2);

% prepare the RIR for blockwise processing
RIRsq = RIR_LR.^2;
RIRsq = [repmat(RIRsq(1,:), N2, 1); RIRsq; repmat(RIRsq(end,:), N2+1, 1)];

% initialize the output and sum of squares
rmsBlock_fs = nan(size(RIR,1), 1);
sq          = sum(RIRsq(1:N,:));

% blockwise estimate
for nn = 1:size(RIRsq,1)-N
    
    rmsBlock_fs(nn,:) = sq;
    sq = sq - RIRsq(nn,:) + RIRsq(nn+N,:);
    
end

rmsBlock_fs = sqrt(rmsBlock_fs / N);

% ---- resample blockwise RMS to the output sample rate HRTF_fs ---- %
if fs ~= HRTF_fs
    t_original = (0:size(RIR_LR,1)-1)' / fs;
    t_target   = (0:1/HRTF_fs:t_original(end))';
    
    rmsBlock = interp1(t_original, rmsBlock_fs, t_target, 'linear');
else
    rmsBlock = rmsBlock_fs;
end

% ---- generate the decay curve ---- %
if strcmpi(tailMethod{1}, 'exact')
    
    % decay curve equals residual RMS energy
    dc = rmsBlock;
    
elseif contains(tailMethod{1}, 'ramp')
    
    % initialze the decay curve
    if numel(tailMethod) == 2
        dc = -inf*ones(numel(rmsBlock),              1);
    else
        dc = -inf*ones(round(tailMethod{3}*HRTF_fs), 1);
    end
    
    % - estimate the apply the fall ramp - %
    % sample of latest included refelctions (plus some breathing time given
    % by the HRTF_delay and half the expected maximum ITD)
    n                  = ceil( max(r.t) * HRTF_fs ) + HRTF_delay + ceil(350e-6*HRTF_fs);
    [pFall, ~, muFall] = polyfit((n:numel(rmsBlock))', db(rmsBlock(n:end)), tailMethod{2});
    fFall              = polyval(pFall, n:numel(dc), [], muFall);
    % apply the fall ramp
    dc(n:end)          = fFall;
    
    % - estimate the rise ramp - %
    % start sample of first sound in HRTF_fs units
    m = ceil( min(r.t)*HRTF_fs ) + HRTF_delay;
    
    if strcmpi(tailMethod{1}, 'ramp_a')
        fRise = linspace(dc(n)-60, dc(n), 10)';
    elseif strcmpi(tailMethod{1}, 'ramp_b')
        fRise = linspace(dc(n)-60, dc(n), n-m+1)';
    elseif strcmpi(tailMethod{1}, 'ramp_c')
        E  = sum(RIR_LR(m:n-1).^2);
        a0 = 10^(dc(n)/20);

        % this should work but seems to be numerically instable
%         if E>0
%             a0   = dc(n);
%             a0   = 10^(a0/20);
%             a0   = a0^2;
%             t    = (m-n) / HRTF_fs;
%             
%             lw    = a0*t * exp( a0*t/(2*E) ) / (2*E);
%             delta = lambertw( lw );
%             delta = a0*t - 2*E * delta;
%             delta = delta / (2*E*t);
%         
%             fRise = a0 * exp( ((m-n):0)'/HRTF_fs * delta );
%             fRise = db(fRise);
%         else
%             fRise = dc(n);
%         end
        
        % this is a bruit force approximation
        if E==0
            fRise = dc(n);
        else
            % initial parameters
            L     = 0;                              % start level of the rise ramp
            t     = (m-n-1:0)'/(m-n-1);             % time axis
            fRise = t .* L + dc(n);                 % rise ramp
            Erise = sum(10.^(fRise(1:end-1) / 10)); % energy of the rise ramp
            
            % increment of the start level
            if Erise > E
                dL = -10e3;
            else
                dL =  10;
            end
            
            % iteratively change the start level until the energy matches
            % the target within .5 dB
            while abs(db(Erise/E)) > .5
                L     = L + dL;
                fRise = t .* L + dc(n);
                Erise = sum(10.^(fRise(1:end-1) / 10));
                
                if Erise == E
                    break
                elseif (E > Erise && dL < 0) || (E < Erise && dL > 0)
                    dL = -.5 * dL;
                end
                
            end
        end

    elseif strcmpi(tailMethod{1}, 'ramp_d')
        if any(abs(rmsBlock(m:n)))
            [pRise, ~, muRise] = polyfit((m:n)', db(rmsBlock(m:n)), tailMethod{2});
            fRise              = polyval(pRise, m:n, [], muRise);
        else
            fRise              = dc(n);
        end
    end
    
    % write rise ramp to final vector
    dc(n-numel(fRise)+1:n) = fRise;
    
    % de-logarithmize
    dc = 10.^(dc/20);
    
else
    error('composeParametricBRIR:Input', ['tailMethod ''' tailMethod{1} ''' not implemented.'])
end

% ---- generate binaurl noise and apply the decay curve ---- %
% white gaussion noise with binaural coherence (assuming perfectly diffuse
% soundfield)
LR = randn(numel(dc), 2);
LR = binauralCoherence(LR, fs, N);

LR = LR .* dc;
    
% check energy in parametric late reverbeartion vs. actual reverberation
Nmismatch  = floor( size(RIR_LR,1)*HRTF_fs/fs );
try
    LRmismatch = db( mean( rms( LR(1:Nmismatch,:) ) ) / rms(RIR_LR) );
catch
    LRmismatch = NaN;
    warning('composeParametricBRIR:Input', 'Duration of late reverberation shorter than RIR. Increase duration if this is unintended (param2 for tailMethods ''ramp'' and ''ramps''.')
end

% resample final  to input sample rate decay curve
if HRTF_fs ~= fs
    t_target = (0:numel(dc)-1)'/HRTF_fs;
    dc_fs    = interp1(t_target, dc, t_original, 'linear', 'extrap'); % 'extrap' catches the case where the duration of the reverberation is shorter than the input, if the user passes a duration. It's not needed otherwise.
else
    dc_fs = dc;
end

%% ------------------------------------------------------------- get output
% ER might be shorter in case 'ramp' or 'ramps' was used and a duration for
% LR.S
ER(end+1:size(LR,1),:) = 0;

% get other params
BRIR       = ER + LR;
LRdecay    = dc;
LRdecayOrg = dc_fs;

%% ------------------------------------------------------------------- plot

if doPlot
    
    tEarly = floor(max(r.t*1e2))*10+50;
    
    dNorm = max(RIR_ER);
    
    ph = newFig(12,24, show_plot);
    ph.Name = 'BRIR';
    subplot(4,1,1)
        plotET(RIR_LR/dNorm, fs, false, [.7 .7 .7], '-')
        plotET(RIR_ER/dNorm, fs, false, 'k',        '-', false, 90)
        YLim = get(gca, 'YLim');
        plotET(rmsBlock/dNorm, HRTF_fs, false, 'o', '-', false, YLim)
        plotET(dc/dNorm,       HRTF_fs, false, 'p', '-', false, YLim)
        legend('Residual SRIR', 'Early reflections', 'Residual RMS energy', 'Late reverb. decay', 'location', 'NorthEast')
        title(['RIR' plotTitle])
    subplot(4,1,2)
        plotET(RIR_LR/dNorm, fs, false, [.7 .7 .7], '-')
        plotET(RIR_ER/dNorm, fs, false, 'k',        '-', false, 60)
        YLim = get(gca, 'YLim');
        plotET(rmsBlock/dNorm, HRTF_fs, false, 'o', '-', false,      YLim)
        plotET(dc/dNorm,       HRTF_fs, false, 'p', '-', [0 tEarly], YLim)
        title(['RIR early' plotTitle])
    subplot(4,1,3)
        plotET(ER(:,1)/dNorm, HRTF_fs, false, 'b',        '-', [0 tEarly], YLim)
        plotET(LR(:,1)/dNorm, HRTF_fs, false, [.7 .7 .7], '-', false,      YLim)
        plotET(ER(:,1)/dNorm, HRTF_fs, false, 'b',        '-', [0 tEarly], YLim)
        title(['Left parametric BRIR early' plotTitle])
        legend('Early reflections', 'Late reverberation', 'location', 'NorthEast')
    subplot(4,1,4)
        plotET(ER(:,2)/dNorm, HRTF_fs, false, 'r',        '-', [0 tEarly], YLim)
        plotET(LR(:,2)/dNorm, HRTF_fs, false, [.7 .7 .7], '-', false,      YLim)
        plotET(ER(:,2)/dNorm, HRTF_fs, false, 'r',        '-', [0 tEarly], YLim)
        title(['Right parametric BRIR early' plotTitle])
        legend('Early reflections', 'Late reverberation', 'location', 'NorthEast')
        
end