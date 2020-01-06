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
function plotET(data, fs, normalize, color, line, xLim, yLim)

if normalize
    data = data / max(abs(data(:)));
end

data             = db(data);
data(data==-Inf) = -300;

t = (0:size(data,1)-1)' / fs * 1e3;
if ~exist('xLim', 'var')
    xLim = [0 t(end)];
end
if islogical(xLim)
    xLim = [0 t(end)];
end

% colors from AKcolors.m from AKtools
if ischar(color)
    colors.rgb    = [0 0 0; 1 1 1; 0 .4470 .7410; .85 .225 .098; .929 .694 .125; .494 .184 .556; .266 .674 .188; .301 .745 .933; .635 .078 .184; .8 .2 .8; .1 .8 .8; 1 .54 0];
    colors.string = 'kwbrypgldmco';
    color         = colors.rgb(strfind(colors.string, color),:);
end

hold on
plot(t, data, line, 'color', color, 'LineWidth', 1)
hold off

box on
grid on
set(gca, 'TickLength', [0 0])
xlim(xLim)
if exist('yLim', 'var')
    if numel(yLim)==1
        tmp_max = max(data(:));
        tmp_max = tmp_max - mod(tmp_max, 10)+10;
        if ceil(max(data(:))) == tmp_max
            tmp_max = tmp_max + 10;
        end
        
        yLim = [tmp_max-yLim tmp_max];
    end
    ylim(yLim)
end

xlabel 'Time in ms'
ylabel 'Amplitude in dB'
title  'Energy Time Curve (ETC)'