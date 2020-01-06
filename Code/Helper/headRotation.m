% Re-licensed from AKtools. See AKroomSimulationRotation.m inside
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
function [azRot, elRot] = headRotation(az, el, rotAz, rotEl)

% check input format
try
    az    = reshape(az,    [numel(az)    1]);
    el    = reshape(el,    [numel(az)    1]);
    rotAz = reshape(rotAz, [1 numel(rotAz)]);
    rotEl = reshape(rotEl, [1 numel(rotAz)]);
catch
    error('headRotation:Input', 'Input data has wrong size: ''az'' and ''el'' must be [N x 1], ''rotAz'' and ''rotEl'' must be [1 x M].')
end

% convert input points to Cartesian coordinates
[x, y, z] = sph2cart(az/180*pi, el/180*pi, ones(size(az)));
X         = [x'; y'; z'];

% apply the inverse rotation
azRot = nan(size(az,1), size(rotAz,2));
elRot = azRot;

for nn = 1:numel(rotAz)
    % rotation matrix about z-Axis
    rAz = [cosd(rotAz(nn)) -sind(rotAz(nn)) 0;...
           sind(rotAz(nn))  cosd(rotAz(nn))  0;...
           0                0                1];

    % rotation matrix about y-Axis
    rEl = [cosd(-rotEl(nn)) 0 sind(-rotEl(nn));...
           0                 1 0;...
           -sind(-rotEl(nn)) 0 cosd(-rotEl(nn))];
    
    % apply inverse rotation to the incoming/outgoing sound paths
    % (rAz^-1 = rAz', and (rAz*rEl)' = rEl'*rAz')
    Xr = (rAz*rEl)' * X;

    % convert to spherical coordinates
    [azRot(:,nn), elRot(:,nn)] = cart2sph(Xr(1,:)', Xr(2,:)', Xr(3,:)');
end

% rad to deg
azRot = azRot/pi*180;
azRot = mod(azRot, 360);
elRot = elRot/pi*180;