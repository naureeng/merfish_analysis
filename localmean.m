function mean = localmean(f, nhood)
% computes the mean at the centre of every neighbourhood in f defined by nhood

if nargin == 1
    nhood = ones(3) / 9;
else
    nhood = nhood / sum(nhood(:));
end
mean = imfilter(f, nhood, 'replicate');