% Runs the parametric encoding/decoding engine proposed in [1]
%
% a.jueterbock@tu-berlin.de, Audio Communication Group TU Berlin &
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin &
% Microsoft Research, Redmond, USA
%
% [1] Tobias Jüterbock, Fabian Brinkmann, Hannes Gamper, 
%     Nikunj Raghuvanshi, and Stefan Weinzierl
%      'Spatio-Temporal Windowing for Encoding Perceptually Salient Early 
%       Reflections in Parametric Spatial Audio Rendering'
%     JAES - https://www.aes.org/e-lib/browse.cfm?elib=22239

%   Copyright 2022 Microsoft Corporation
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
close all; clear; clc

% REQUIREMENTS:
% DSP System Toolbox
% Signal Processing Toolbox 
% AKtools (https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/)

% add paths
addpath('Code/')
addpath('Code/Helper/')
addpath SRIRs

% check for HRIR
if ~exist('FABIAN_HRIR_measured_HATO_0', 'file')
    getHRIR
end

%%
% --------------------------------------------------------- select the room
room.volume        = 'small';   % 'small', 'medium', or 'large'
room.reverberation = 'dry';     % 'dry', 'medium', or 'wet'
room.position      = 'corner';  % 'center', 'wall', or 'corner'
room.name          = ['ROOM-'           room.volume        ...
                      ' REVERBERATION-' room.reverberation ...
                      ' POSITION-'      room.position];


% ------------------------------------ parameters for detecting reflections
%                                        (see detectReflections.m for help)
    setup.timeMax     = 85;
    setup.timeRange   = [.5 .8];
    setup.angleRange  = [1 1 1 1 1];
    setup.thEcho_lat  = [-0.06   0 1  0 .5 0 0]; 
    setup.thEcho_pol  = [-0.06   0 1  0 .5 0 0]; 
    setup.thMask_lat  = [-1 -17.25 1 10 .5 0 1]; 
    setup.thMask_pol  = [-1 -17.25 1 10 .5 0 1]; 
    setup.eEcho       = [1 .05];
    setup.eMask       = [1 .35];


% ----------------------- parameters for reducing the number of reflections
%                                        (see reduceReflections.m for help)
setup.reduceN       = 4;
setup.reduceMethod  = 'loudest';


% -------------------------------- parameters for parametric BRIR rendering
%                                    (see composeParametricBRIR.m for help)
setup.diffuseN      = 256;
setup.diffuseMethod = {'ramp_c' 1};

% sampling rate of FABIAN HRTFs (do not change)
setup.HRTF_fs       = 44100;

% --------------------------------------------------------------- load data
data = load(room.name);

%% -----------------------------------------------------check for toolboxes
if ~license('test','signal_toolbox')
    warning("It seems that the Signal Processing Toolbox is not installed.")
end
if exist('hor2sph','file') ~= 2
    warning("It seems that AKtools is not inside the Matlab search path. It is available from https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/")
end

%% ----------------------------------------------------- detect reflections
fprintf('--------- %s ---------\n', room.name)
fprintf('Detect reflections\n')
[er.r, er.t, er.isAudible, er.tEcho, er.tMask] = detectReflections(data.doa, data.rir, data.fs, setup.timeMax, setup.timeRange, setup.angleRange, setup.thEcho_lat, setup.thEcho_pol, setup.thMask_lat, setup.thMask_pol, setup.eEcho, setup.eMask, data.t_mix, true);

% ------------------------------------------------------ reduce refelctions
fprintf('Reduce reflections\n')
[rr.r, rr.id, rr.isAudible] = reduceReflections(er.r, setup.reduceN, setup.reduceMethod, data.fs, room.name, true);

% -------------------------------------------------- render parametric BRIR
fprintf('Render parametric BRIR\n')
rng(7)
[brir.h, brir.er, brir.lr, brir.LRmismatch] = composeParametricBRIR(data.rir, data.doa, data.fs, rr.r, rr.isAudible, [0 0], setup.diffuseN, setup.diffuseMethod, false, false, room.name, true);


fprintf('done\n\n')