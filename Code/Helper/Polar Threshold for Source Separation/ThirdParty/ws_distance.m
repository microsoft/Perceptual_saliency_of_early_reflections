function wsd = ws_distance(u_samples, v_samples, p)
% WS_DISTANCE 1- and 2- Wasserstein distance between two discrete 
% probability measures 
%   
%   wsd = WS_DISTANCE(u_samples, v_samples) returns the 1-Wasserstein 
%   distance between the discrete probability measures u and v 
%   corresponding to the sample vectors u_samples and v_samples
%
%   wsd = WS_DISTANCE(u_samples, v_samples, p) returns the p-Wasserstein 
%   distance between the discrete probability measures u and v
%   corresponding to the sample vectors u_samples and v_samples. 
%   p must be 1 or 2.
%
% from https://github.com/nklb/wasserstein-distance

if ~exist('p', 'var')
    p = 1;
end

u_samples_sorted = sort(u_samples(:));
v_samples_sorted = sort(v_samples(:));

if p == 1
    
    all_samples = unique([u_samples_sorted; v_samples_sorted], 'sorted');
    
    u_cdf = find_interval(u_samples_sorted, all_samples(1:end-1)) ...
        / numel(u_samples);
    v_cdf = find_interval(v_samples_sorted, all_samples(1:end-1)) ...
        / numel(v_samples);
    
    wsd = sum(abs(u_cdf - v_cdf) .* diff(all_samples));
    
elseif p == 2
    
    u_N = numel(u_samples);
    v_N = numel(v_samples);    
    all_prob = unique([(0:u_N) / u_N, (0:v_N) / v_N], 'sorted').';
    
    u_icdf = u_samples_sorted(fix(all_prob(1:end-1) * u_N) + 1);
    v_icdf = v_samples_sorted(fix(all_prob(1:end-1) * v_N) + 1);
    
    wsd = sqrt(sum((u_icdf-v_icdf).^2 .* diff(all_prob)));
    
else
    
    error('Only p=1 or p=2 allowed.')
    
end
end

function idx = find_interval(bounds, vals)
% Given the two sorted arrays bounds and vals, the function 
% idx = FIND_INTERVAL(bounds, vals) identifies for each vals(i) the index 
% idx(i) s.t. bounds(idx(i)) <= vals(i) < bounds(idx(i) + 1).

m = 0;
bounds = [bounds(:); inf];
idx = zeros(numel(vals), 1);

for i = 1:numel(vals)
    while bounds(m+1) <= vals(i)
        m = m + 1;
    end
    idx(i) = m;
end
end
% MIT License
% 
% Copyright (c) 2020 Niklas Kolbe
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
