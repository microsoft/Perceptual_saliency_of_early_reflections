% plotSRIR(RIR, DOA, plane, fs, t_max, pr, pn, type, ls, c, cbar)
% plots the spatial room impulse response (SRIR) in the current axes or 
% creates a figure if none exists. DIR is plotted on the lateral or polar
% angle vs. time.
%
% I N P U T
% RIR   - [N 1] vector of the room impulse response. Will be normalized to
%         0 dB for plotting.
% DOA   - x/y/z vector pointing to the direction of arrival. Size [N 3]
%         where N is the length in samples
% plane - 'lat' for plotting the lateral angle (default)
%         'pol' for plotting the polar angle
%         'rir' draws only the RIR with the color indicating the level and
%               using type = 'line'.
%         'all' generates a 3 by one subplot layout showing the RIR and the
%               SRIR in the lateral and polar plane
% fs    - sampling frequency in Hz (default = 44100)
% t_max - maximum time in s up to which the data is plotted. Pass 'range'
%         to plot everything according to the plot range, which is the
%         default (see 'pr' below).
% pr    - plot range in dB. Sets the range of the z-axis (default = 30).
% pn    - true   normalize the plot to 0 dB (default).
%         false  no normalization.
% type  - 'scatter'  draws a marker at the location of each sample and
%                    indicates the level by the color and line height (3D)
%                    (default)
%         'line'     draws a line plus marker at the location of each
%                    sample and indicates the level by the color and line
%                    height (3D, slower than scatter)
% ls    - a three element cell array {marker linewidth markersize} that
%         gives the linespecifications. The default {'.' 1 10} plots lines
%         with a width of 2 and a dot marker with a size of 10
%         (type 'doc linespec' for more info).
% c     - rbg vector of size [1 3] giving the plot color or string
%         specifying a colormap ('R', 'G', or 'B'). The resolution
%         of the colormap is 1 dB. (default = 'G' for Green)
% cbar  - true to show colorbar. false to hide colorbar (default = true)
%
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
function plotSRIR(RIR, DOA, plane, fs, t_max, pr, pn, type, ls, c, cbar)

% default values
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('t_max', 'var')
    t_max = 'range';
end
if ~exist('pr', 'var')
    pr = 30;
end
if ~exist('pn', 'var')
    pn = true;
end
if ~exist('plane', 'var')
    plane = 'lat';
end
if ~exist('type', 'var')
    type = 'scatter';
end
if ~exist('ls', 'var')
    ls = {'.' 1 10};
end
if ~exist('c', 'var')
    c = 'Greens';
end
if ~exist('cbar', 'var')
    cbar = true;
end

pr = abs(pr);

% plot all
if strcmpi(plane, 'all')
    % make figure if none exists
    if isempty(findall(0,'Type','Figure'))
        newFig(12,18)
    end

    subplot(3,1,1)
        plotSRIR(RIR, DOA, 'rir', fs, t_max, pr, pn, type, ls, c, cbar)
    subplot(3,1,2)
        plotSRIR(RIR, DOA, 'lat', fs, t_max, pr, pn, type, ls, c, cbar)
    subplot(3,1,3)
        plotSRIR(RIR, DOA, 'pol', fs, t_max, pr, pn, type, ls, c, cbar)
    return
end

% set range of the RIR
h           = db(abs(RIR(:,1)));
h(isinf(h)) = nan;
if pn
    h = h - max(h);
end
h(h<-pr+max(h)) = NaN;

% truncate to t_max
if ischar(t_max)
    n_max = find( ~isnan(h), 1, 'last' );
    t_max = n_max / fs;
    t_max = t_max - rem(t_max, 10e-3) + 10e-3;
else
    n_max = min( round(t_max*fs)+1, size(h,1) );
end
DOA = DOA(1:n_max,:);
h   = h(1:n_max);

% get normalized version to pick the color
hn  = h - max(h);


% get polar angle from DOA: interpret horizontal polar format as rotated 
% spherical coordinates with negative azimuth direction (taken from SOFA
% Matlab API, sph2hor.m)
[pol,nlat] = cart2sph(DOA(:,1),DOA(:,3),-DOA(:,2));
pol = rad2deg(pol);
lat = rad2deg(-nlat);
% adjust polar angle range
pol = mod(pol+90,360)-90;

if strcmpi(plane, 'rir')
    ylab = 'Amplitude in dB';
elseif strcmpi(plane, 'lat')
    p_ang  = lat;
    ylab   = 'Lateral angle in degree';
    yli    = [-90 90];
    ytick  = -90:22.5:90;
    ytilab = {-90 '' -45 '' 0 '' 45 '' 90};
elseif strcmpi(plane, 'pol')
    p_ang = pol;
    ylab  = 'Polar angle in degree';
    yli   = [-90 270];
    ytick = -90:45:270;
    ytilab = {-90 '' 0 '' 90 '' 180 '' 270};
else
    error('plotSRIR:input', 'plane must be ''lat'', ''pol'', or ,''rir''.')
end

% generate the colormap
d2 = round(pr/2);
d2 = d2 + mod(d2,2);
if isnumeric(c)
    c = repmat(reshape(c, [1 3]), pr+1, 1);
else
    % intial colormap (variation from ColorBrewer)
    if strcmpi(c, 'b')
        c = [247,251,255;222,235,247;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107.01] / 255;
    elseif strcmpi(c, 'r')
        c = [255,245,240;254,224,210;252,187,161;252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13.01] / 255;
    elseif strcmpi(c, 'g')
        c = [247,252,245;229,245,224;199,233,192;161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27.01] / 255;
    end
    
    % interpolate colormap
    c = [interp1(1:9, c(:,1)', linspace(1,9,round(pr+d2)), 'linear')' ...
         interp1(1:9, c(:,2)', linspace(1,9,round(pr+d2)), 'linear')' ...
         interp1(1:9, c(:,3)', linspace(1,9,round(pr+d2)), 'linear')'];
    
    % take desired sections
    c = c(d2/2:pr+d2/2,:);
end

% plot
% make figure if none exists
if isempty(findall(0,'Type','Figure'))
    newFig(12,6)
end
hold on

if strcmpi(plane, 'rir')
    hp = h;
    hp(isnan(hp)) = ceil(max(h))-pr-.5;
    hp(end)       = nan;
    patch((0:n_max-1)/fs*1e3, hp, hp, 'EdgeColor', 'interp', 'Marker', 'none', 'LineWidth', ls{2})
    colormap(c)
elseif strcmpi(type, 'line')
    for nn = 1:numel(lat)
        if ~isnan(h(nn))
            plot3([nn-1 nn-1]/fs*1e3, [p_ang(nn) p_ang(nn)], [-pr+max(h)-2 h(nn)], 'color', c(round(hn(nn)+pr+1),:), 'LineWidth', ls{2})
            plot3((nn-1)/fs*1e3, p_ang(nn), h(nn), ls{1}, 'color', c(round(hn(nn)+pr+1),:), 'LineWidth', ls{2}, 'MarkerSize', ls{3})
        end
    end
elseif strcmpi(type, 'scatter')
    for nn = 1:numel(p_ang)
        if ~isnan(h(nn))
            plot3((nn-1)/fs*1e3, p_ang(nn), h(nn), ls{1}, 'color', c(round(hn(nn)+pr+1),:), 'LineWidth', ls{2}, 'MarkerSize', ls{3})
        end
    end
else
    error('plotSRIR:input', 'type must be ''line'' or ''scatter''.')
end

% format axis
hold off
xlim([0 t_max*1e3])
xlabel 'Time in ms'

if strcmpi(plane, 'rir')
    ylabel(ylab)
    ylim([ceil(max(h))-pr ceil(max(h))])
    
    set(gca, 'YTick', flip(0:-10:-pr)+round(max(h)))
else
    ylabel(ylab)
    zlabel 'Amplitude in dB'
    
    ylim(yli)
    zlim([ceil(max(h))-pr ceil(max(h))+1])
    set(gca, ...
        'YTick', ytick, 'YTickLabel', ytilab, ...
        'ZTick', flip(0:-10:-pr)+round(max(h)))
    
    view([0 90])
end

box on
grid on
rotate3d on

if cbar
    colormap(c)
    cb = colorbar;
    if strcmpi(plane, 'rir')
        set(cb, 'Ticks', flip(0:-10:-pr)+ceil(max(h)), 'Limits', [ceil(max(h))-pr ceil(max(h))])
    else
        set(cb, 'Ticks', flip( (pr:-10:0)/pr ), 'TickLabels', flip(0:-10:-pr)+round(max(h)))
    end
    ylabel(cb, 'Amplitude in dB')
end