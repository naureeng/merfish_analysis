function I = percentile2i(h, P)
% function that given a percentile I and a histogram h, computes an
% intensity I representing the Pth percentile and returns the value in I

if P < 0 || P > 1
    error('The percentile must be in range [0,1].');
end

% normalise h
h = h/sum(h);

% cumulative distribution
C = cumsum(h);

idx = find(C >= P, 1, 'first');

% subtract 1 from idx because indeixing starts at 1 but intensities start
% at 0, plus normalise range to [0 1]
I = (idx - 1) / (numel(h) - 1);

end