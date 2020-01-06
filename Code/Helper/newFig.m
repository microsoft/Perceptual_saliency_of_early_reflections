% Creates an empty figure. Re-licensed from AKtools. See AKf.m inside
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
function hFigureHandle = newFig(fWidth, fHeight, num)


if nargin == 2
    
    hFigureHandle = figure;
    
elseif nargin == 3
    
    if isnumeric(num)
        hFigureHandle = figure(num);
    elseif num
        hFigureHandle = figure;
    else
        hFigureHandle = figure('visible','off');
    end
    
end

% get screensize in centimeter
pixpercm = get(0,'ScreenPixelsPerInch')/2.54;
screen = get(0,'ScreenSize')/pixpercm;
% get position for figure
left   = max([(screen(3)-fWidth)/2 0]);
bottom = max([(screen(4)-fHeight)/2 0]);

set(hFigureHandle,'PaperUnits', 'centimeters');
set(hFigureHandle,'Units', 'centimeters');

% paper size for printing
set(hFigureHandle, 'PaperSize', [fWidth fHeight]);
% location on printed paper
set(hFigureHandle,'PaperPosition', [.1 .1 fWidth-.1 fHeight-.1]);
% location and size on screen
set(hFigureHandle,'Position', [left bottom fWidth fHeight]);

% set color
set(hFigureHandle, 'color', [1 1 1])

if ~nargout
    clearvars hFigureHandle
end
