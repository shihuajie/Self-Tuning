%% create_multires
%% Input
%   img                 the image source
%   multires            the number of multi-resolution levels
%   patch_size          the size of a basic patch
%   multires_weight     the weight of the multires levels
%% Output
%   out_row             the output row of multiresolution levels
function [ out_row ] = create_multires( img, multires, patch_size, multires_weight )
    out_row = cell(1, 1 + multires);
    out_row{1} = img;
    P = patch_size;
    % P(r) = patch_size * 3^(r-1) = 3 * P(r-1), P(1) = patch_size
    % TL(r) = -patch_size * (3^r-1)/2
    if nargin < 4
        multires_weight = 1;
    end
    R = numel(multires_weight);
    for r = 1:multires
        h = fspecial('average', [P, P]);
        if R > 1
            w = multires_weight(r);
        else
            w = multires_weight;
        end
        out_row{1 + r} = w * imfilter(img, h, 'replicate', 'same');
        P = P * 3;
    end
end

