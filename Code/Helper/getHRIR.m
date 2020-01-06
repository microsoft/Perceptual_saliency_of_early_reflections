% Downloads and extracts FABIAN HRIRs
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
function getHRIR

if ~exist('FABIAN_HRIR_measured_HATO_0.mat', 'file')
    
    download = input('The FABIAN HRIRs are missing. Download 100 MB automatically? (Y/N): ', 's');
    
    if strcmpi(download, 'y')
        disp('Downloading https://depositonce.tu-berlin.de/bitstream/11303/6153.3/8/FABIAN_HRTFs_NeutralHeadOrientation.zip')
        websave(fullfile('Code', 'Helper', 'HRIR.zip'), 'https://depositonce.tu-berlin.de/bitstream/11303/6153.3/8/FABIAN_HRTFs_NeutralHeadOrientation.zip')
        
        disp('Unzipping')
        unzip(fullfile('Code', 'Helper', 'HRIR.zip'), fullfile('Code', 'Helper'))
        
        disp('Cleaning data')
        movefile(fullfile('Code', 'Helper', '0 HRIRs neutral head orientation', 'SphericalHarmonics', 'FABIAN_HRIR_measured_HATO_0.mat'), fullfile('Code', 'Helper', 'FABIAN_HRIR_measured_HATO_0.mat'))
        rmdir(fullfile('Code', 'Helper', '0 HRIRs neutral head orientation'), 's')
        delete(fullfile('Code', 'Helper', 'HRIR.zip'))
        
        fprintf('\n\n')
    else
        error('Download and unzip https://depositonce.tu-berlin.de/bitstream/11303/6153.3/8/FABIAN_HRTFs_NeutralHeadOrientation.zip inside your Matlab search path.')
    end
end