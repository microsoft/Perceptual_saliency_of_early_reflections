% Re-licensed from AKtools. See AKisht.m inside
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
function y = isht(f_nm, doIFFT, multi, SHTmode, is_even, compact, SHmode)

if ~exist('SHTmode', 'var')
    SHTmode = 'db_unwrap';
end
if ~exist('is_even', 'var')
    is_even = 1;
end
if ~exist('compact', 'var')
    if size(f_nm,3) == 1 && any(strcmpi(SHTmode, {'db_unwrap' 'abs_unwrap' 'real_imag'}))
        compact = true;
    else
        compact = false;
    end
end
if ~exist('SHmode', 'var')
    SHmode = 'complex';
end

% ---------------------------- get Ynm matrix if not passed to the function
if size(multi, 2)==2
    N   = sqrt(size(f_nm,1))-1;
    Ynm = sh(N, [], multi(:,1), multi(:,2), SHmode);
else
    Ynm = multi;
end

% ---------------------------------------------------- inverse SH transform
Y = zeros(size(Ynm,1), size(f_nm,2), size(f_nm,3));

for n = 1:size(f_nm,2)
    for c = 1:size(f_nm,3)
        Y(:,n,c) = Ynm * f_nm(:,n,c);
    end
end

% transpose (could have been done above, but this makes it more
% illustrative)
tmp = Y;
Y   = zeros(size(f_nm,2), size(Ynm,1), size(f_nm,3));
for c = 1:size(f_nm,3)
    Y(:,:,c) = tmp(:,:,c).';
end

clear tmp

% -------------------------------- post process x to match desired SHT mode
if compact && any(strcmpi(SHTmode, {'db_unwrap' 'abs_unwrap' 'real_imag'}))
    Y   = cat(3, real(Y), imag(Y));
end

switch SHTmode
    case 'db_unwrap'
        y = 10.^(Y(:,:,1)/20) .* exp(1j*Y(:,:,2));
    case 'abs_unwrap'
        y = Y(:,:,1) .* exp(1j*Y(:,:,2));
    case 'complex'
        y = Y;
    case 'real_imag'
        y = Y(:,:,1) + 1j*Y(:,:,2);
    case 'db'
        y = 10.^(Y/20);
    case 'abs'
        y = Y;
    otherwise
        error(['SHTmode ''' SHTmode ''' unknown'])
end

% ------------------------------------ transfer into time domain if desired
if doIFFT
    y = ifft(single2bothSidedSpectrum(y, is_even), 'symmetric');
end
