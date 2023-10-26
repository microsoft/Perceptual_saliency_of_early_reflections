% [th, th_pol] = polarThreshold(pol_target, lat_ref, pol_ref, width, ...
%                                        depth, transition, doPlot)
% calculates a threshold for the audibility of concurrent sources in
% dependency of the polar angle (assuming an interaural polar coordinate
% system).
%

% TODO: change to lp_target, lp_ref because they're stored together anyway

%
% I N P U T
% pol_target - polar angles of the targets (e.g. concurrent sources,
%              early reflections). Scalar, vector, or [] in which case only
%              th_pol is calculated (see output).
% lat_ref    - lateral angle of the reference (e.g. direct sound)
% pol_ref    - polar angle of the reference (e.g. direct sound)
% width      - scalar that determines the width of the range where high
%              masking is assumed.
%              1: according to values extracted from
%              calculatePolarThreshold.m
%             <1: smaller range
%             >1: larger range
%              (default = 1)
% depth      - difference in dB between low and high masking (default = 10)
% transition - amount of transition between between low and high masking. 
%              0: no transition
%              1: widest transition
%              (default = 0.5)
% doPlot     - 1 to plot polar Threshold curve of all targets 
%              (default = 0)
%
% O U T P U T
% th         - treshold curve at the angles given in pol_target of size
%              [N x 1], where N is the number of target angles.
% th_pol     - two element vector containing the polar limits of the
%              masking threshold
% th_width   - 
% 
% Written by Tobias Jüterbock, Audio Communication Group, TU Berlin
% Structure adapted from lateralThreshold.m, written by
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin &
% Microsoft Research, Redmond, USA
%
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
% 
function [th, th_pol,th_width] = polarThreshold(pol_target, lat_ref, pol_ref, width, ...
                                                depth, transition, doPlot)
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
if ~exist('doPlot', 'var')
    doPlot = false;
end
% check input
depth = abs(depth);
% lat_ref = abs(lat_ref);

if numel(pol_ref) > 1 || numel(transition) > 1 || numel(depth) > 1
    error('concurrentSourceAudibility:Input', ...
          'pol_ref, transition, and depth must be scalars not vectors.')
end

if pol_ref < -90 || pol_ref > 270 || any(pol_target<-90) || any(pol_target>270)
    error('concurrentSourceAudibility:Input', ...
          'pol_target and pol_ref must be in the range -90 <= pol <= 270')
end

if transition > 1 || transition < 0
    transition = max(transition, 0);
    transition = min(transition, 1);
    warning('concurrentSourceAudibility:Input', ...
            'transition was clipped to 0<=transition<=1.')
end


%% Load Polar threshold for source separation
% spatial masking thresholds for detecting two simultaneous sources
% 1st dimension: polar angle of first sound
% 2nd dimension: absolute lateral angle of first sound
% 3rd dimension: 1: polar threshold for positive polar opening angle
%                2: polat threshold for negative polar opening angle
TH = load('Helper/polarThreshold_full_all_VPs_absoluteLatAngle_lat90added_a396ae7cf8afe4c00387625fa804ebe0.mat',...
          'pol_ref','lat_ref_abs','polarThreshold');
th_ref = TH.polarThreshold;
th_pol_ref = TH.pol_ref;
th_lat_ref_abs = TH.lat_ref_abs;

%% interpolate width at reference position
% convert meshgrid to Nx2 grid (N = number of grid points):
[g_lat,g_pol] = meshgrid(th_lat_ref_abs,th_pol_ref);
g_lat = reshape(g_lat,[],1);
g_pol = reshape(g_pol,[],1);
% convert to azimuth-elevation:
[g_az, g_el] = hor2sph(g_lat,g_pol);
% convert reference position to azimuth and elevation
[az_ref, el_ref] = hor2sph(abs(lat_ref),pol_ref);

% interpolate polar threshold width with spherical spline interpolation
% 1: positive width
% 2: negative width
% tic;
% th_width = [
%             AKsphSplineInterp(g_az,g_el,width* ...
%              reshape(th_ref(:,:,1),[],1),az_ref,el_ref, 1, 0.001, 'deg', 0)
%             AKsphSplineInterp(g_az,g_el,width* ...
%              reshape(th_ref(:,:,2),[],1),az_ref,el_ref, 1, 0.001, 'deg', 0)
%             ];
% toc;tic
% TODO: this is linear interpolation, values differ by a few degrees but
% it is 20 times faster than spherical spline interpolation
th_width = [
            griddata(g_az, g_el, width*reshape(th_ref(:,:,1),[],1), ...
                     az_ref, el_ref,'linear'); 
            griddata(g_az, g_el, width*reshape(th_ref(:,:,2),[],1), ...
                     az_ref, el_ref,'linear');
            ];
% toc;
        
% clip values to -180 < th < 180:
th_width(1) = min(th_width(1), 179.9999);
th_width(2) = max(th_width(2),-179.9999);

% set the characteristic points of the masking/threshold curve
% th_lat = lat_ref;
th_pol = pol_ref + th_width;
% wrap to the interval -90<=pol<270
th_pol = mod(th_pol+90,360)-90;

clear az_ref el_ref g_az g_el g_lat g_pol

%% Calculate Masking threshold curve
if ~isempty(pol_target) && depth ~= 0
    % get the masking/threshold curve
    if transition == 0
        % ---- binary masking threshold ---- %
        th = zeros(numel(pol_target), 1);
        if th_pol(1) >= th_pol(2) 
            %find the targets between the reference position and
            %threshold width
            id = pol_target <= th_pol(1) & pol_target >= th_pol(2);  
            % set the threshold of all other targets to -depth
            th(~id) = -depth;
        elseif th_pol(1) < th_pol(2) 
            % find the targets between the reference position and
            % threshold width
            id = pol_target <= th_pol(1) | pol_target >= th_pol(2); 
            % set the threshold of all other targets to -depth
            th(~id) = -depth;
        elseif th_pol(1) == th_pol(2)
            if th_pol(1) == pol_ref && th_pol(2) == pol_ref
                print('polar threshold 0')
                id = pol_target == pol_ref;
                th(~id) = -depth;

            elseif th_width(1) == 180 && th_width(2) == -180
                print('polar threshold 180')
            else
                print('check values, thresholds are weird.')
            end
        end
    else
        % ----- tanh transitioned masking threshold ----- %
        % express difference of pol_target to pol_ref and the th_pol in units
        % of th_width
        d_pol     = nan(numel(pol_target), 1).';
        
        %unwrap pol_target to the interval pol_ref-180 <= pol_target < pol_ref+180:
        pol_target_unwrapped = pol_target;
        id = pol_target < pol_ref - 180;
        pol_target_unwrapped(id) = pol_target_unwrapped(id) + 360;
        id = pol_target >= pol_ref + 180;
        pol_target_unwrapped(id) = pol_target_unwrapped(id) - 360;
        
        %unwrap th_pol to the interval pol_ref-180 <= th_pol < pol_ref+180:
        th_pol_unwrapped = th_pol;
        id = th_pol < pol_ref - 180;
        th_pol_unwrapped(id) = th_pol(id) + 360;
        id = th_pol >= pol_ref + 180;
        th_pol_unwrapped(id) = th_pol(id) - 360;
            
        % difference of pol_target to pol ref in units of th_width:
        id = pol_target_unwrapped >= pol_ref;
        d_pol(id) = (th_pol_unwrapped(1) - pol_target_unwrapped(id) ) ./ th_width(1);
        id = pol_target_unwrapped < pol_ref;
        d_pol(id) = (th_pol_unwrapped(2) - pol_target_unwrapped(id) ) ./ th_width(2);
       
        % adjust transition so it decreases for large polar thresholds and
        % is 0 for a polar threshold of 180 degrees
        a = 160;
        b = 180;
        if th_width(1) > a && abs(th_width(2)) > a
            x = max(abs(th_width));
            transition = transition * ...
             (1 + ((x-b) .* heaviside(x-b) - (x-a) .* heaviside(x-a)) / (b-a));
        end
        
        % adjust depth so it decreases for large polar thresholds and
        % is 0 for a polar threshold of 180 degrees
        a = 160;
        b = 180;
        if th_width(1) > a && abs(th_width(2)) > a
            x = min(abs(th_width));
            depth = depth * ...
             (1 + ((x-b) .* heaviside(x-b) - (x-a) .* heaviside(x-a)) / (b-a));
        end
            
        % get threshold curve by interpolating between 0 dB and -depth dB using
        % a tangens hyperbolicus function
        th = tanh(1.5*d_pol*pi/transition) * depth/2 - depth/2;
    end
    % Plot
    if doPlot
    plot(pol_target,th);
    xlim([min(pol_target),max(pol_target)])
    ylim([-abs(depth)-0.1*abs(depth),0.1*abs(depth)])
    xline(pol_ref);
    xlabel('polar target angle $[^\circ]$','interpreter','latex')
    ylabel('polar threshold')
    end
elseif ~isempty(pol_target) && depth == 0
    th = zeros(numel(pol_target), 1);
else
    th = [];
end
end

