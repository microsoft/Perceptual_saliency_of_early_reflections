% Runs the parametric encoding/decoding engine proposed in [1]
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin &
% Microsoft Research, Redmond, USA
%
% [1] Fabian Brinkmann, Hannes Gamper, Nikunj Raghuvanshi, and Ivan Tashev
%     'Towards encoding perceptually salient early reflections for
%      parametric spatial audio rendering.' 148th AES Convention (Accepted)

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
close all; clear; clc
    
% add paths
addpath(genpath('Code'))
addpath SRIRs

% check for HRIR
if ~exist('FABIAN_HRIR_measured_HATO_0', 'file')
    getHRIR
end

% --------------------------------------------------------- select the room
room.volume        = 'small';   % 'small', 'medium', or 'large'
room.reverberation = 'dry';     % 'dry', 'medium', or 'wet'
room.name          = ['ROOM-' room.volume ' REVERBERATION-' room.reverberation];


% ------------------------------------ parameters for detecting refelctions
%                                        (see detectReflections.m for help)
setup.timeMax       = 85;
setup.timeRange     = [.5 1];
setup.angleRange    = [1 6 1];
setup.thEcho        = [-0.06   0 1  0 .5 0 0];
setup.thMask        = [-1.00 -10 1 10 .5 0 1];
setup.eEcho         = [1 .05];
setup.eMask         = [1 .35];


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


%% ----------------------------------------------------- detect reflections
fprintf('--------- %s ---------\n', room.name)
fprintf('Detect reflections\n')
[er.r, er.t, er.isAudible, er.tEcho, er.tMask] = detectReflections(data.doa, data.rir, data.fs, setup.timeMax, setup.timeRange, setup.angleRange, setup.thEcho, setup.thMask, setup.eEcho, setup.eMask, data.t_mix, room.name);

% ------------------------------------------------------ reduce refelctions
fprintf('Reduce reflections\n')
[rr.r, rr.id, rr.isAudible] = reduceReflections(er.r, setup.reduceN, setup.reduceMethod, data.fs, room.name, true);

% -------------------------------------------------- render parametric BRIR
fprintf('Render parametric BRIR\n')
rng(7)
[brir.h, brir.er, brir.lr, brir.LRmismatch] = composeParametricBRIR(data.rir, data.doa, data.fs, rr.r, rr.isAudible, [0 0], setup.diffuseN, setup.diffuseMethod, false, false, room.name, true);


fprintf('done\n\n')