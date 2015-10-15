function [im_out,im_mask] = fix_downsampled_mask( input_im, input_mask)
%input_mask(input_mask < 0.01) = 0;

TT = 0.9;
ch_num = size(input_im,3);
%input_mask(input_mask > TT) = 1;
mask_3ch = repmat(input_mask, [1 1 ch_num]);

im_out = input_im ./ ( 1 - mask_3ch );
im_out(mask_3ch == 1) = 0;
im_mask = input_mask;
im_mask(im_mask < TT) = 0;
im_mask(im_mask ~= 0) = 1;
%im_out(repmat(input_mask, [1 1 ch_num]) == 1) = 0;
