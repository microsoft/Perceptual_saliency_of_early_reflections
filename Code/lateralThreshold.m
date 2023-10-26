% [th, th_lat] = lateralThreshold(lat_target, lat_ref, width, depth, transition)
% claculates a threshold for the audibility of concurrent sources in
% dependency of the lateral angle (assuming an interaural polar coordinate
% system).
%
% I N P U T
% lat_target - lateral angles of the targets (e.g. concurrent sources,
%              early reflections). Scalar, vector, or [] in which case only
%              th_lat is calculated (see output).
% lat_ref    - lateral angle of the reference (e.g. direct sound)
% width      - scalar that determines the width of the range where high
%              masking is assumed.
%              1: according to values extracted from [1, Fig. 3(e)]
%             <1: smaller range
%             >1: larger range
%              (default = 1)
% depth      - difference in dB between low and high masking (default = 10)
% transition - amount of transition between between low and high masking. 
%              0: no transition
%              1: widest transition
%              (default = 0.5)
%
% O U T P U T
% th         - treshold curve at the angles given in lat_target of size
%              [N x 1], where N is the number of target angles.
% th_lat     - two element vector containing the lateral limits of the
%              masking threshold
%
% [1] V. Best, A. van Schaik, and S. Carlile (2004): 'Separation of
%     concurrent broadband sound sources by human listeners.' J. Acoust.
%     Soc. Am. 115(1):324-336.
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
function [th, th_lat] = lateralThreshold(lat_target, lat_ref, width, depth, transition)

% default parameters
if ~exist('width', 'var')
    width = 1;
end
if ~exist('depth', 'var')
    depth = 10;
end
if ~exist('transition', 'var')
    transition = .5;
end

% check input
depth = abs(depth);

if numel(lat_ref) > 1 || numel(transition) > 1 || numel(depth) > 1
    error('concurrentSourceAudibility:Input', 'lat_ref, transition, and depth must be scalars not vectors.')
end

if abs(lat_ref) > 90 || any(abs(lat_target)>90)
    error('concurrentSourceAudibility:Input', 'lat_target and lat_ref must be in the range -90 <= lat <= 90')
end

if transition > 1 || transition < 0
    transition = max(transition, 0);
    transition = min(transition, 1);
    warning('concurrentSourceAudibility:Input', 'transition was clipped to 0<=transition<=1.')
end

% spatial masking thresholds for detecting two simultaneous sources
% 1 column: absolute lateral angle of first sound
% 2 column: offset towards median plane
% 3 column: offset away from median plane
% Values in columns 2 and 3 were manually extracted from [1, Fig. 3(e)] and
% denote the point close to 100% detection. The value in the third column
% at 67.5 deg was linearly interpolateed from it's neighbors.
th_ref = [
                 0    -15 15
                 22.5 -20 30
                 45   -30 35
                 67.5 -38 40
                 90   -45 45
                ];

% interpolate width of the spatial masking threshold to reference position
th_width = [
              interp1( th_ref(:,1), th_ref(:,2)*width, abs(lat_ref) ) ...
              interp1( th_ref(:,1), th_ref(:,3)*width, abs(lat_ref) )
            ];

% set the characteristic points of the masking/threshold curve
if lat_ref >= 0
    th_lat = lat_ref + [th_width(2) th_width(1)];
else
    th_lat = lat_ref - th_width;
end

if ~isempty(lat_target) && depth ~= 0
    % get the masking/threshold curve
    if transition == 0
        % ---- binary masking threshold ---- %
        th = zeros(numel(lat_target), 1);
        id = lat_target <= th_lat(1) & lat_target >= th_lat(2);
        th(~id) = -depth;
    else
        % ----- tanh transitioned masking threshold ----- %
        % express difference of lat_target to lat_ref and the th_lat in units
        % of th_width
        d_lat     = nan(numel(lat_target), 1);
        
        id        = lat_target >= lat_ref;
        if lat_ref > 0
            d_lat(id) = ( th_lat(1) - lat_target(id) ) ./      th_width(2);
        else
            d_lat(id) = ( th_lat(1) - lat_target(id) ) ./ abs( th_width(1) );
        end
        
        id        = lat_target <  lat_ref;
        if lat_ref > 0
            d_lat(id) = ( lat_target(id) - th_lat(2) ) ./      th_width(2);
        else
            d_lat(id) = ( lat_target(id) - th_lat(2) ) ./ abs( th_width(1) );
        end
        % get threshold curve by interpolating between 0 dB and -depth dB using
        % a tanges hyperbolicus function
        th = tanh(1.5*d_lat*pi/transition) * depth/2 - depth/2;
    end
elseif ~isempty(lat_target) && depth == 0
    th = zeros(numel(lat_target), 1);
else
    th = [];
end