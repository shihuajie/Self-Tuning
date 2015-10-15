function [out_im, out_mask] = random_merge(im1, im2, blend_alpha, patch_size)
nCh = size(im1,3);
szX = size(blend_alpha,1);
szY = size(blend_alpha,2);
out_im = zeros(szX, szY, nCh);
dimX = ceil(szX / patch_size);
dimY = ceil(szY / patch_size);
rnd = rand(dimX, dimY);

posX = rand(dimX, dimY);
posY = rand(dimX, dimY);
blend_im_rsz = imresize(blend_alpha,size(rnd),'nearest');


for i = 1 : dimX
    dstXS = (i - 1) * patch_size + 1;
    dstXE = min(i * patch_size, szX);
    lX = dstXE - dstXS;
    for j = 1 : dimY
        
        dstYS = (j - 1) * patch_size + 1;
        dstYE = min(j * patch_size, szY);
        lY = dstYE - dstYS;
        if rnd(i,j) < blend_im_rsz(i,j)
            rnd(i,j) = 0;
            srcXS = floor((size(im1,1) - patch_size) * posX(i,j)) + 1;
            srcYS = floor((size(im1,2) - patch_size) * posY(i,j)) + 1;
            out_im(dstXS : dstXE, dstYS : dstYE, :) = im1(srcXS: srcXS + lX, srcYS: srcYS+lY,:);
        else
            rnd(i,j) = 1;
            srcXS = floor((size(im2,1) - patch_size) * posX(i,j)) + 1;
            srcYS = floor((size(im2,2) - patch_size) * posY(i,j)) + 1;
            out_im(dstXS : dstXE, dstYS : dstYE, :) = im2(srcXS: srcXS + lX, srcYS: srcYS+lY,:);
        end
    end
end
out_mask = imresize(rnd, [szX, szY],'nearest');
