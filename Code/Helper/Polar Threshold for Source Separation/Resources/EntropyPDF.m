function E = Entropy(p)
%Entropy Calculates the Entropy of each column of p
%
% Calculated according to "A Mathematical Theory of Communication"
% By C. E. Shannon 
%
% I N P U T  : 
% p          : Array, discrete probability density functions in each column
%              (each column has to sum to 1 and contain values between 0 and 1)
%
% O U T P U T:
% E          : Entropy of each column
% 
% Written by Tobias Jüterbock

E = -sum(p .*log2(p), 1,'omitnan');
end