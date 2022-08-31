close all
clear all

%% Load polar threshold

filepath = './_Results/Threshold';
prefix   = 'polarThreshold_all_VPs';
checksum = 'a396ae7cf8afe4c00387625fa804ebe0';
filename = [prefix '_' checksum];
load (fullfile(filepath,filename));

%% Average positive and negative lateral angles
idx = lat_ref >= 0;
idxLat0 = find(idx,1);

lat_ref_abs = lat_ref(idxLat0:end);
assert(lat_ref_abs(1) == 0, ...
    'This routine expects lat angles to be symmetric around lat=0')

polarThreshold_pos = 0.5 * (polarThreshold(:,idxLat0:end ,1) + ...
                            polarThreshold(:,idxLat0:-1:1,1));
polarThreshold_neg = 0.5 * (polarThreshold(:,idxLat0:end ,2) + ...
                            polarThreshold(:,idxLat0:-1:1,2));

polarThreshold = nan(numel(pol_ref),numel(lat_ref_abs),2);
polarThreshold(:,:,1) = polarThreshold_pos;
polarThreshold(:,:,2) = polarThreshold_neg;

%% add polarThreshold = +-180 for lateral angle 90:
idx = numel(lat_ref_abs);
lat_ref_abs(idx+1) = 90;
polarThreshold(:,idx+1,1) =  180 * ones(numel(pol_ref),1);
polarThreshold(:,idx+1,2) = -180 * ones(numel(pol_ref),1);

%% save to file
filename = [prefix '_absoluteLatAngle_lat90added_' checksum '.mat'];

save(fullfile(filepath,filename),'lat_ref_abs','pol_ref','polarThreshold');