function [ res ] = labelize( img, N, N0 )
%LABELIZE Generate a label map
    if nargin < 2 || N < 2
        N = 64;
    end
    if nargin < 3
        N0 = 1;
    end

    m = min(img(:));
    M = max(img(:));
    res = int32((img - m) / (M - m) * (N - 1)) + N0;
end

