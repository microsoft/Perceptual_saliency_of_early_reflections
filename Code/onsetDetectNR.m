% ons = onsetDetectNR(h, fs, fLP, noiseLevel)
% estimates the onset of the first sound arrival in impulse responses
% according to [1], Eq. (15).
%
% I N P U T
% h          - impulse response of size [N C] with N being the number of
%              samples and C the number of channels
% fs         - sampling rate in Hz (default = 44100)
% fLP        - cut-off frequency (-3 dB) of first order Butterworth lowpass
%              that is used to smooth the cumulative energy
%              (default = fs/4)
% noiseLevel - estimated level in dB of the noise relative to the absolute
%              maximum value of h (default = -60)
%
% O U T P U T
% ons        - estimated onset in fractional samples
%
%
% [1] N. Raghuvanshi and J. Snyder (2018): "Parametric directional coding
%     for precomputed sound propagation." ACM Trans. Graph. 37(4), Article
%     108.
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
function ons = onsetDetectNR(h, fs, fLP, noiseLevel)

% ------------------------------------------------------ set default values
if ~exist('fs', 'var')
    fs = 44100;
end

if ~exist('fLP', 'var')
    fLP = fs/4;
end

if ~exist('noiseLevel', 'var')
    noiseLevel = -60;
end

% make sure it works with any use input
noiseLevel = -abs(noiseLevel);

% -------------------------------------------------------- onset estimation

% acount for two samples that get lost due to the ratio and derivate in
% [1], Eq. (15) - see below (not mentioned in [1])
h = [zeros(2, size(h,2)); h];

% --- normalize h to the absolute maximum across channels to make sense of
% the noiseLevel, which is assumed to be a peak2tail value
h = h / max(abs(h(:)));

% --- calucalte E(t) according to [1], Eq. (15)
E_t = cumsum(h.^2);
if fLP
    % lowpass
    Hd = fdesign.lowpass('N,F3dB', 1, fLP/fs*2);
    Hd = design(Hd, 'butter');
    if ~isstable(Hd)
        error('onsetDetectNR:Input', 'Filter instable. Try to increase fLP.')
    end
    E_t = filter(Hd, E_t);
end

% ---- calculate D(t) according to [1], Eq. (15)
% first derivative considering the noise level
D_t        = E_t(2:end,:) ./ (E_t(1:end-1,:) + 10^(noiseLevel/10));
% remove everything below the noiseLevel (not described in [1], but might
% make it more robust in some cases)
D_t(D_t<1) = 1;
% second derivative to make D(t) peak at sudden energy jumps
D_t        = D_t(2:end,:) -  D_t(1:end-1,:);
% controll peakedness as in [1]
D_t        = D_t.^2;

% calucalte tau_0 based on D(t)
t     = (0:size(h,1)-3)';
tau_0 = cumsum(t.*D_t) ./ cumsum(D_t);

% last value gives the estimated onset
% NOTE: This is implemented slightly different in [1] due to realtime
% demands. But as mentioned in [1] this is the way to go if not doing
% real-time
ons = tau_0(end,:);