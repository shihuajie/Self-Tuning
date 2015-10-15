function [] = imagesc_bar( I, aspect, cbar )
%IMAGESC_BAR Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 3
        cbar = 1;
    end
    if nargin < 2
        aspect = 1;
    end

    % show stuff
    imagesc(I);
    if aspect
        axis image;
    end
    if cbar
        colorbar;
    end

end

