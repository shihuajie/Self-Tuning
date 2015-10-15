%%
% randomInit - initialize a target image by shuffling random patches of a
% source target
%
%% Input
%   trg_im      the target image
%   src_im      the source texture
%   P_x         the X block size
%   P_y         the Y block size (P_x)
%   int_loc     whether to use integer location for the blocks
%
%% Output
%   out_im      an image of the target size, with randomly
%               shuffled patches from the source texture
function [ out_im ] = random_init( trg_im, src_im, P_x, P_y, int_loc )

    if nargin < 4
        P_y = P_x;
        int_loc = 0;
    elseif nargin < 5
        int_loc = 0;
    end

    [trgSzX, trgSzY, nCh] = size(trg_im);
    out_im = zeros(trgSzX, trgSzY, nCh);
    [srcSzX, srcSzY, ~] = size(src_im);
    
    % grid of blocks
    trgDimX = ceil(trgSzX / P_x);
    trgDimY = ceil(trgSzY / P_y);

    % blocks position
    if int_loc % integer locations
        srcDimX = floor(srcSzX / P_x);
        srcDimY = floor(srcSzY / P_y);
        posX = randi(srcDimX - 1, trgDimX, trgDimY);
        posY = randi(srcDimY - 1, trgDimX, trgDimY);
    else
        posX = rand(trgDimX, trgDimY);
        posY = rand(trgDimX, trgDimY);
    end

    % for each block of the target
    for i = 1 : trgDimX
        % beginning and end of the block (X)
        dstXS = (i - 1) * P_x + 1;
        dstXE = min(i * P_x, trgSzX);
        lX = dstXE - dstXS; % effective size of the block (X)
        
        for j = 1 : trgDimY
            % beginning and end of the block (Y)
            dstYS = (j - 1) * P_y + 1;
            dstYE = min(j * P_y, trgSzY);
            lY = dstYE - dstYS; % effective size of the block (Y)
            
            % beginning and end of the block in the source
            if int_loc
                srcXS = posX(i, j) * P_x;
                srcYS = posY(i, j) * P_y;
            else
                srcXS = floor((srcSzX - P_x) * posX(i,j)) + 1;
                srcYS = floor((srcSzY - P_y) * posY(i,j)) + 1;
            end
            out_im(dstXS : dstXE, dstYS : dstYE, :) = src_im(srcXS : srcXS+lX, srcYS : srcYS+lY,:);
        end
    end

end