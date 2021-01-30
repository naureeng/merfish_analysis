function g = localthresh(f, nhood, a, b, meantype)
% thresholds image f by computing local threshold at centre (x,y) of every
% neighbourhood in f.
% nhood defines size of neighbourhood
% 
% segmented image is given by
%       1 if (f > a*sig) and (f > b*mean)
%   g = 
%       0 otherwise

% % initialise
% f = tofloat(f);

% compute local STDs
SIG = stdfilt(f, nhood);
% compute MEAN
if nargin == 5 && strcmp(meantype,'global')
    MEAN = mean2(f);
else
    MEAN = localmean(f, nhood);
end

% obtain segmented image
g = (f > a*SIG) & (f > b*MEAN);

end
