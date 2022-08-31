close all; clear; clc
% requires Parallel Computing Toolbox. If not available, change the parfor
% loop into a for loop

% define grid:
lat_ref = -85:5:85;
pol_ref = -90:5:269;
pol_mu = 4; % stepsize of the opening angle

% HUTUBS database HRTF datasets:
VPs = [1:87, 89:95];

printProgress = false;
savedata = true;

% unique identifier
checksum = DataHash([lat_ref, pol_ref,pol_mu,VPs]);

%% Allocate space

p = zeros(360,360/pol_mu,numel(pol_ref),numel(lat_ref), 'single');

fprintf('Shape of results:\n\n'); disp(size(p));
fprintf('(pol x pol_ref x pol_delta x lat_ref)\n')
assert(360/pol_mu == round(360/pol_mu))

%% Model perceived location of 2 concurrent sources with the baumgartner model

fprintf('Total number of iterations: %i\n', numel(lat_ref))
% parallel loop for each lateral angle:
parfor ll = 1:numel(lat_ref)
    %check for results that have been calculated already:
    lat = lat_ref(ll);
    
    pathname = './_Results/baumgartner/partial_results_for_each_lat/'
    filename = sprintf('baumgartner_all_VPs_lat%i_%s.mat',lat,checksum);
    if isfile(fullfile(pathname,filename))
        % if result exist, load results:
        fprintf('Loading results for lateral angle %i°. \n',lat)
        partialdata = load(fullfile(pathname,filename));
        p(:,:,:,ll) = partialdata.data;
        fprintf('Done loading results for lateral angle %i°. \n',lat)
    else
        % if results don't exist, calculate them:
        tic
        fprintf('Calculating lateral angle %i°. \n',lat)

        for VP = VPs
        tic
        fprintf('Calculating VP %i\n',VP)
        data = calculateBaumgartner(lat, pol_ref, pol_mu, VP, printProgress, false);
        p(:,:,:,ll) = p(:,:,:,ll) + data.p;
        toc
        end

        if savedata
            pathname = './_Results/baumgartner/partial_results_for_each_lat';
            filename = sprintf('baumgartner_all_VPs_lat%i_%s', lat,checksum);
            if ~isfolder(pathname)
                mkdir(pathname)
            end
            parsave(fullfile(pathname,filename), p(:,:,:,ll))
        end
        fprintf('Done calculating lateral angle %i°. \n',lat)
        toc
    end

end

% Average results across subjects
p = p ./ numel(VPs);

%% Write results to file
baumgartner_all_VPs.pol = (-90:269).';
baumgartner_all_VPs.pol_mu = pol_mu;
baumgartner_all_VPs.pol_delta = 0:abs(pol_mu):359;
baumgartner_all_VPs.pol_ref = pol_ref;
baumgartner_all_VPs.lat_ref = lat_ref;
baumgartner_all_VPs.model = 'baumgartner';
baumgartner_all_VPs.VPs = VPs;
baumgartner_all_VPs.p = p;


if savedata
    pathname = './_Results/baumgartner/';
    filename = sprintf('baumgartner_all_VPs_%s',checksum);
    if ~isfolder(pathname)
        mkdir(pathname)
    end
    save(fullfile(pathname,filename), 'baumgartner_all_VPs')
end

