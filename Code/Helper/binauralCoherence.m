% Applies binaural coherence. Re-licensed from AKtools. See
% AKbinauralCoherence inside www.ak.tu-berlin.de/AKtools for a complete and
% documented version.
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
function y = binauralCoherence(x, fs, N)

if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('N', 'var')
    N = 2^12;
end

% in case of non-stationary filtering for generation of the reverb, we want
% to avoid having a bin at fs/2 in our fft. This bin had to be real and
% would commonly be set to 0 (-inf dB) and thus we would have a dip in the
% spectrum...
if ~mod(N, 2)
    N = N+1;
end

% frequency vector
f = (0:fs/N:fs/2)';

% interaural coherence from Borss2009 eq.(8)
gamma    = sin(pi*f/550)./(pi*f/550) .* max([zeros(size(f)) 1-f/2700], [], 2);
gamma(1) = 1;

% filter for coherence adjustment from Borss2009 eq.(17)(18)
H_b          = sqrt(1/2 * (1-sqrt(1-gamma.^2)));
H_b(gamma<0) = -H_b(gamma<0);

H_a = sqrt(1-H_b.^2);

% get double sided spectra and impulse responses
H_b = single2bothSidedSpectrum(H_b, 1-mod(N,2));
H_a = single2bothSidedSpectrum(H_a, 1-mod(N,2));

h_b = ifft(H_b, 'symmetric');
h_a = ifft(H_a, 'symmetric');

h_b = circshift(h_b, [round(N/2), 0]);
h_a = circshift(h_a, [round(N/2), 0]);

% filter
y(:,1) = fftfilt(h_a, [x(:,1); zeros(N-1,1)]) + fftfilt(h_b, [x(:,2); zeros(N-1,1)]);
y(:,2) = fftfilt(h_a, [x(:,2); zeros(N-1,1)]) + fftfilt(h_b, [x(:,1); zeros(N-1,1)]);

% circshift to account for the delay of the filter
y = circshift(y, [-round(N/2) 0]);
y = y(1:end-N+1,:);

