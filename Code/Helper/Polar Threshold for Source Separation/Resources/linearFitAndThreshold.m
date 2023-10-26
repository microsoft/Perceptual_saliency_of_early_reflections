function [out, fitparams] = linearFitAndThreshold(xvalues, yvalues, threshold, regInterval, limit)
% linearFitAndThreshold perform linear fit and threshold the fit
% Written by Tobias Jüterbock
% 
% Performs a linear fit on the x-value interval specified in regInterval
% and calculates the intersection of the fit and the threshold value.
% Limits the output to the interval specified in limit.
%
% I N P U T  :
% xvalues    : vector containing the x values 
% yvalues    : vector containing the y values
% threshold  : threshold value
% regInterval: 2 entry vector containing the regression min and max values
% limit      : limit of the output
%
% O U T P U T:
% out        : intersection of the fit with the threshold value
% fitparams  : fit parameters
%
% Example:
% xvalues = 0:360;
% yvalues = xvalues.^2;
% [out, f] = linearFitAndThreshold(xvalues,yvalues,1E4,[0 100],180)

% fitparams = cell(size(regvalues));
xvalues = xvalues(:);
yvalues = yvalues(:);


% set the x-Interval where linear regression should be calculated
valid = xvalues <= max(regInterval) & xvalues >= min(regInterval);

% Do a linear fit to the data:
% (Limit the quadratic term so the parabolic term is 0)
fitparams = fit(xvalues(valid), yvalues(valid), 'poly1');

% calculate intersection with the threshold value:
out = (threshold-fitparams.p2)/fitparams.p1;

% limit output:
if limit >= 0 
    out = min(out,limit);
    % in case the fit has a slope that would yield to a negative
    % intersection, set out to limit:
    if out < 0 
        out = limit;
    end
elseif limit < 0
    out = max(out,limit);
    % in case the fit has a slope that would yield to a negative
    % intersection, set out to limit:
    if out > 0 
        out = limit;
    end
end
    
end