function [ F ] = smart_features( I0, F0, N, Fproc, Dproc )
%SMART_FEATURES Extract features
%
% SYNOPSIS
%   [ F ] = smart_features( I, D, N )
%
% INPUT
%   - I0    the Lab image
%   - F0    the edge map
%   - N     the number of classes to look for
%   - Fproc the processing applied to each feature channel
%
% OUTPUT
%   - F     the features
%
    
    % simplified image using slic
    I = improc(I0, 'slic:1:20');

    if nargin < 3
        N = 2;
    end
    
    if nargin < 4 || isempty(Fproc)
        Fproc = 'rev,sbwdist,n:0:100';
    end
    if nargin < 5 || isempty(Dproc)
        Dproc = Fproc;
    end
    
    % workspace
    [H, W, ~] = size(I);

    % compute edge mask
    D = imdilate(F0, ones(5), 'same'); % edge mask
    D0 = reshape(D, H * W, 1);

    % color data
    C = reshape(I(:, :, 1:3), H * W, 3);
    % isolate the edge part
    C(repmat(D0, [1, 3]) == 1) = 1000;
    
    % use mean-shift to do clustering
    [cluster_idx, cluster_center] = kmeans(C, N+1, 'distance', 'sqEuclidean', 'Replicates', N+1);
    pixel_labels = reshape(cluster_idx, H, W);
    % imshow(pixel_labels, []); title('Segmentation');
    
    F = zeros(H, W, N+1);
    Dch = 1;
    maxInter = 0;
    for f = 1:N+1
        Fn = double(pixel_labels == f);
        inter = numel(find(Fn == D));
        if inter > maxInter
            maxInter = inter;
            Dch = f;
        end
        F(:, :, f) = Fn; % imdilate(Fn, ones(5), 'same');
        % figure; imshow(F(:, :, f)); title(sprintf('Feature %d', f));
    end
    
    % remove edge mask from the non-edge channels, intersect it with edge
    % channel
    for f = 1:N+1
        Fn = F(:, :, f);
        if f == Dch
            F(:, :, f) = improc(double(Fn & D), Dproc);
        else
            F(:, :, f) = improc(double(Fn & ~D), Fproc);
        end
        
    end
    1;
    
    % remove edge feature
%     if N > 1
%         F = F(:, :, [1:Dch-1, Dch+1:size(F, 3)]);
%     end
end

