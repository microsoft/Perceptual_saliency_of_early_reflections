close all; clear; clc

% add project files
addpath(genpath('Resources'))

% add anechoic audio, FABIAN HRTFs, generated SRIRs, ...
addpath(genpath(fullfile('..', 'Data')))

% add single scripts in third party foler
addpath ThirdParty

% add AKtools (Fabian Brinkmann)
cd(fullfile('ThirdParty', 'AKtools'))
AKtoolsStart
cd(fullfile('..', '..'))

% add SDM toolbox (Sakarai Tervo)
addpath(fullfile('ThirdParty', 'SDM Toolbox'))
addpath(fullfile('ThirdParty', 'SDM Toolbox', 'Data'))

% add AMToolbox
addpath(genpath(fullfile('ThirdParty', 'amtoolbox-full-0.10.0')))

% add SOFiA toolbox
addpath(genpath(fullfile('ThirdParty', 'SOFiA-master')))

% add ITA toolbox
addpath(genpath(fullfile('ThirdParty', 'ITAtoolbox')))