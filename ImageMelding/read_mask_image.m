%reads the input image and wherever the color is magenta it assings it to be a hole!
function [im_out im_mask invalid_mask] = read_mask_image(f_name, in_scale,thr)
im_out = im2double(imread(f_name));
imR = im_out(:,:,1);

my_ones = ones(size(imR));
im_mask_R = zeros(size(imR));
inds = find(imR >= thr);
im_mask_R( inds) = my_ones(inds);
im_mask_RN = zeros(size(imR));
im_mask_RN(imR <= 1 - thr) = 1;
imG = im_out(:,:,2);
im_mask_G = zeros(size(imG));
inds = find(imG <=  1 - thr);
im_mask_G( inds) = my_ones(inds);
im_mask_GN = zeros(size(imG));
im_mask_GN(imG >= thr) = 1;
imB = im_out(:,:,3);
im_mask_B = zeros(size(imB));
inds = find(imB >= thr);
im_mask_B( inds) = my_ones(inds);


im_mask = double(im_mask_R & im_mask_B & im_mask_G);
invalid_mask = double(im_mask_RN & im_mask_B & im_mask_GN);
inds = find(repmat(im_mask, [1 1 3]));
invalid_inds = find(repmat(invalid_mask, [1 1 3]));
im_out(invalid_inds) = 1000;
im_out(inds) = 0;
im_out = imresize(im_out,in_scale);
im_mask = imresize(im_mask,in_scale);
[im_out im_mask] = fix_downsampled_mask(im_out,im_mask);
my_ones = ones(size(im_mask));
inds = find(im_mask ~= 0);
im_mask(inds) = my_ones(inds);
