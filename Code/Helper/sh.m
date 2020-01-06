% Re-licensed from AKtools. See AKsh.m inside
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
function [Ynm, N, M] = sh(n, m, az, el, mode)

% --- input parsing ---
if numel(n)~=1
    error('sh:Input', 'n must be a scalar value')
elseif rem(n,1)
    error('sh:Input', 'n must be an integer value')
end

if~exist('mode', 'var')
    mode = 'complex';
end

az = reshape(az, [numel(az) 1]);
el = reshape(el, [numel(el) 1]);

if size(az,1)~=size(el,1) || size(az,2)~=size(el,2)
    error('sh:Input', 'az and el must be of same size')
end

% coordinate conversion
az  = az/180*pi;
el  = el/180*pi;

% --- get combinations of n and m to compute ---    
% number of spherical harmonics, N,M,az matrices
[N, M, Nsh] = getNM(n);

N  = repmat(N', [numel(az) 1]);
M  = repmat(M', [numel(az) 1]);
az = repmat(az, [1 Nsh]);


clear nn num_sh pos_start

% --- calculate spherical harmonics ---
% scalar (dependend on n,m; independend from az, el)

% check if all factorials can be computed
if any( [N-M N+M] > 170 )
    error('sh:input', 'SH basis functions can only be computed up to order 85. Otherwise the factorials can not easiliy be computed with built in variable types.')
end

a = sqrt((2*N+1)./(4*pi) .* factorial(N-M)     ./factorial(N+M));

% azimuthal changing part of spehrical harmonics
e = exp(1j*M.*az);

% elevation dependend part of spherical harmonics
if Nsh > 1
    l = ones(size(N));
    % get combinations
    for nn = 1:n
        pos_start = (nn)^2+1;   % number of previously calculated n,m-combinations
        
        % calculate legendre functions for positive m
        l(:,pos_start+nn:pos_start+2*nn) = legendre(nn, cos(el))';
        if strcmpi(mode, 'complex')
            % legendre functions for negative m from [2], eq. (6.31)
            mm = 1:nn;
            lScale = (-1).^mm .* factorial(nn-mm)./factorial(nn+mm);
            l(:,pos_start:pos_start+nn-1) = fliplr(l(:,pos_start+nn+1:pos_start+2*nn) .* repmat(lScale, [size(az,1) 1]));
        else
            l(:,pos_start:pos_start+nn-1) = fliplr( l(:,pos_start+nn+1:pos_start+2*nn) );
        end
    end
else
    l = legendre(n, cos(el))';
    if ~isempty(m)
        l = l(:,abs(m)+1);
        if m < 0 && strcmpi(mode, 'complex')
            % legendre function for negative m from [2], eq. (6.31)
            l = l * (-1).^-m .* factorial(n+m)./factorial(n-m);
        end
    end
end

% get spherical harmonics
Ynm = a .* l .* e;

% get corresponding n,m-vectors
N = N(1,:);
M = M(1,:);

end

function [N, M, Nsh] = getNM(n)

% number of spherical harmonics, N,M,az matrices
Nsh = (n+1)^2;
N   = zeros(Nsh, 1);
M   = N;

% get combinations
for nn = 1:n
    pos_start = (nn)^2+1;   % number of previously calculated n,m-combinations
    N(pos_start:pos_start+2*nn) = nn;
    M(pos_start:pos_start+2*nn) = -nn:nn;
end
end