close all; clear; clc

% add project files
addpath(genpath('Resources'))

% add single scripts in third party foler
addpath ThirdParty

%% -----------------------------------------------------check for toolboxes
if ~license('test','signal_toolbox')
    warning("It seems that the Signal Processing Toolbox is not installed.")
end
if ~license('test','distrib_computing_toolbox')
    warning("It seems that the Parallel Computing Toolbox is not installed.")
end
if ~license('test','curve_fitting_toolbox')
    warning("It seems that the Curve Fitting Toolbox is not installed.")
end
if exist('hor2sph','file') ~= 2
    warning("It seems that AKtools is not inside the Matlab search path. It is available from https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/")
end
if exist('sofia_fliege','file') ~= 2
    warning("It seems that the SOFiA sound field analysis toolbox is not inside the Matlab search path. It is available from https://audiogroup.web.th-koeln.de/SOFiA_wiki/WELCOME.html")
end
if exist('AKhrirDatabase', 'file') ~= 2
    warning("It seems that the HUTUBS HRTF database is not inside the Matlab search path. It is available from https://dx.doi.org/10.14279/depositonce-8487")
end
if exist('baumgartner2014', 'file') ~= 2
    warning("It seems that the Auditory Modeling Toolbox is not inside the Matlab search path. It is available from https://amtoolbox.org/")
end