function threshold = thresholdLateralIncrease(threshold_median, lat)
%thresholdLateralIncrease threshold increase for nonzero lateral angles
% Written by Tobias Jüterbock
%
% INPUT :
% threshold_median: threshold in the median plane
% lat             : lateral angle(s)
% 
% OUTPUT :
% threshold       : corrected threshold for each lateral angle
% 
% for constant polar angle differences, the elevation difference decreases
% with increasing lateral angle. 
%
% 
% lat = 0:89;
% pol = 1;
% [~, el]      = hor2sph(lat, pol);
% plot(lat,1./el)
% hold on
% plot(lat,1./(cos(rad(lat))))

threshold = threshold_median ./ cos(rad(lat));

end

