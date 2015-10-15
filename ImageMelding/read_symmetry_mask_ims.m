function [im1 im1_mask im1_src im2_src blend_alpha hole_mask src2_mask] = read_symmetry_mask_ims(mix_im_f, hole_im_f, src2_hole_f, src1_f, src2_f)
im1_src = im2double(imread(src1_f));
im2_src = im2double(imread(src2_f)); 
hole_mask = rgb2gray(imread(hole_im_f));
src2_mask = rgb2gray(imread(src2_hole_f));
[im1 im1_mask] = read_mask_image(mix_im_f,1,0.8);

% x2 = bwdist(src2_mask);
% x2(im1_mask == 0) = 0;
% xg = imfilter(hole_mask,[0 -1 0;-1 4 -1; 0 -1 0]);
% l = xg > 0;
% myMin = min(x2(l));
% x2(x2 > myMin) = myMin;
% x2 = x2 / myMin;
% x2(hole_mask == 0) = 1;
% 
% blend_alpha = x2;

m1 = bwdist(src2_mask);
m2 = bwdist(1 - hole_mask);
blend_alpha = m1 ./ (m1 + m2);

hole_mask(hole_mask ~= 0) = 1;
src2_mask(src2_mask ~= 0) = 1;
hole_mask = double(hole_mask);
src2_mask =  1 - double(src2_mask);