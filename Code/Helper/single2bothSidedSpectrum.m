% Re-licensed from AKtools. See AKsingle2bothSidedSpectrum.m inside
% www.ak.tu-berlin.de/AKtools for a complete and documented version.
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
function both_sided = single2bothSidedSpectrum(single_sided, is_even)

if ~exist('is_even', 'var')
    is_even = 1;
end

if is_even>1
    is_even = 1-mod(is_even,2);
end

N = size(single_sided,1);

if is_even
    % make sure that the bin at nyquist frequency is real
    % (there might be rounding errors that produce small immaginary parts)
    single_sided(end,:,:) = real(single_sided(end,:,:));
    % mirror the spectrum
    both_sided = [single_sided; flipud(conj(single_sided(2:N-1,:, :)))];    
else
    % mirror the spectrum
    both_sided = [single_sided; flipud(conj(single_sided(2:N,:, :)))];
end

% make sure that the bin at 0 Hz is real
% (there might be rounding errors that produce small immaginary parts)
both_sided(1,:,:) = real(both_sided(1,:,:));