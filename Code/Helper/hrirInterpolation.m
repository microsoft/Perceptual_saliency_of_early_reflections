% Re-licensed from AKtools. See AKhrirInterpolation.m inside
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
function [l, r, az, el, HATO] = hrirInterpolation(az, el, HATO, type)

%% ------------------------------------------------------ 1. HRIR parameter

% resolution of available HATOs
dHATO = 10;
% range of available HATOs
rHATO = [310 50];
% number of frequencies coded by SH coefficients
Nfreq = 129;
% Number of SH coefficients per frequency
Nnm = (35+1)^2;
% number of HRIR samples
Nsamp = 256;
% number of grid points in impulse response data
Ngrid = 11950;


%% --------------------------------------------------------- 2. input check
% separate data type and interpolation domain
type = strsplit(type, '_');
for nn = 1:numel(type)
    if any( strcmpi( type{nn}, {'measured' 'modeled'} ) )
        type_acquisition = type{nn};
    elseif any( strcmpi( type{nn}, {'sh' 'ir'} ) )
        type_domain = type{nn};
    elseif any( strcmpi( type{nn}, {'hrir' 'dir'} ) )
        type_data = type{nn};
    else
        error('hrirInterpolation:type', [type{nn} 'is not a valid flag for ''type'''])
    end
end

if ~exist('type_acquisition', 'var')
    type_acquisition = 'measured';
end
if ~exist('type_domain', 'var')
    type_domain = 'sh';
end
if ~exist('type_data', 'var')
    type_data = 'hrir';
end

clear type

% number of outpt HRIRs
Nhrir = max([numel(el) numel(az) numel(HATO)]);

clear id


%% ---------------------------------------------- 3. Format input arguments
az   = reshape(mod(az,360), numel(az), 1);
el   = reshape(el, numel(el), 1);
HATO = reshape(mod(HATO,360), numel(HATO), 1);

if numel(az) ~= Nhrir
    az = repmat(az(1), Nhrir, 1);
end
if numel(el) ~= Nhrir
    el = repmat(el(1), Nhrir, 1);
end
if numel(HATO) ~= Nhrir
    HATO = repmat(HATO(1), Nhrir, 1);
end

%% -------------------- 4. List source positions and HATO for interpolation
% clip to range of available HATOs
HATO(HATO>rHATO(2) & HATO <=180) = 50;
HATO(HATO>180 & HATO <rHATO(1))  = 310;

% get HATOs for interpolation (eq. (4) in [1])
HATOint      = zeros(Nhrir, 2);
HATOint(:,1) = mod(floor(HATO/dHATO)*dHATO, 360);
HATOint(:,2) = mod((floor(HATO/dHATO)+1)*dHATO, 360);

% don't interpolate HATO if data exists
HATOint(~mod(HATO,dHATO),2) = NaN;

% get needed source positions (eq. (8) in [1])
azInt = mod(HATOint-repmat(HATO, [1 2])+repmat(az, [1 2]), 360);

% get HATOs needed for interpolation
HATOunique = unique(HATOint(:));
HATOunique = HATOunique(~isnan(HATOunique));


%% ----------------------------------------------------------- 5. Load data   
% needed HATOs
hato = HATOunique;

% allocate memory
data_l = zeros(Nnm, Nfreq, numel(hato));
data_r = data_l;

% load
for nn = 1:numel(hato)
    try
        load(['FABIAN_HRIR_' type_acquisition '_HATO_' num2str(hato(nn))], 'SH');
        data_l(:,:,nn) = SH.coeffLeft;
        data_r(:,:,nn) = SH.coeffRight;
    catch
        error('toolbox:dependencies', 'This part of needs the FABIAN HRIR data set which is available from: https://dx.doi.org/10.14279/depositonce-5718.2 inside the Matlab search path.');
    end
end

%% ----------------------------------------- 6. interpolate source position
% allocate space for HRTFs needed for interpolation
lInt = zeros(Nfreq, Nhrir, 2);
rInt = lInt;

for nn = 1:numel(hato)
    for ii = 1:2
        id = HATOint(:,ii)==hato(nn);
        if any(id)
            lInt(:,id,ii) = isht(data_l(:,:,nn), false, [azInt(id,ii) 90-el(id)], 'complex');
            rInt(:,id,ii) = isht(data_r(:,:,nn), false, [azInt(id,ii) 90-el(id)], 'complex');
        end
    end
end

clear tmp
%% ---------------------------------------------------- 7. interpolate HATO
% only interpolate when necessary
idInt = ~isnan(HATOint(:,2));

% allocate space for final HRIRs
l = zeros(Nfreq, numel(HATO));
r = l;

if any(idInt)
    % get weights for interpolation (eq. (5), and (4) in [1])
    HATOweight = acosd(cosd(HATOint(idInt,:)-repmat(HATO(idInt), [1 2])));
    HATOweight = HATOweight ./ repmat(sum(HATOweight,2), [1 2]);

    HATOweightA = repmat(HATOweight(:,2)', [Nfreq 1]);
    HATOweightB = repmat(HATOweight(:,1)', [Nfreq 1]);


    % interpolate HATO separately for magnitude and unwrapped phase
    l(:,idInt) = ( abs(lInt(:,idInt,1)).*HATOweightA + abs(lInt(:,idInt,2)).*HATOweightB ) .* ...
                   exp(1j * ( unwrap(angle(lInt(:,idInt,1))).*HATOweightA + unwrap(angle(lInt(:,idInt,2))).*HATOweightB ) );

    r(:,idInt) = ( abs(rInt(:,idInt,1)).*HATOweightA + abs(rInt(:,idInt,2)).*HATOweightB ) .* ...
                   exp(1j * ( unwrap(angle(rInt(:,idInt,1))).*HATOweightA + unwrap(angle(rInt(:,idInt,2))).*HATOweightB ) );
end

% don't interpolare HATO if not neccessary           
l(:,~idInt) = lInt(:,~idInt,1);
r(:,~idInt) = rInt(:,~idInt,1);

%% ----------------------------------------------------------- 8. get HRIRs
% set bin at fs/2 to 0
l(end,:) = 0;
r(end,:) = 0;

% get HRIRs
l = ifft(single2bothSidedSpectrum(l), 'symmetric');
r = ifft(single2bothSidedSpectrum(r), 'symmetric');

% get DIRs
if strcmpi(type_data, 'dir')
    load FABIAN_CTF %#ok<LOAD>
    
    l = fftfilt(ctf, l);
    r = fftfilt(ctf, r);
end

% WARNING: an earlier version of AKtools reversed the phase of the HRIRs
% because AKsht and AKist used the ' operator for transposing a complex
% vector. However, Matlab calculates the complex conjugate transpose in
% this case. The current version now uses the .' operator.
if strcmpi(type_domain, 'sh')
    if sum(abs(l(1:Nsamp/2,1))) < sum(abs(l(Nsamp/2:end,1)))
        l = flipud(l);
        r = flipud(r);
    end
end