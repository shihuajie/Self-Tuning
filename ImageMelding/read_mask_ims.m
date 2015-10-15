function [im1 im1_mask im1_src im2_src blend_alpha hole_mask src2_mask] = read_mask_ims(im_src_f, im_trg_f, src_msk_f, synth_msk_f)
im1_src = im2double(imread(im_src_f));
im2_src = im2double(imread(im_trg_f)); 
src2_mask = double(rgb2gray(imread(src_msk_f)));
synth_mask = double(rgb2gray(imread(synth_msk_f)));
im1_mask = double(synth_mask);

hole_mask = src2_mask + synth_mask;
hole_mask(hole_mask ~= 0) = 1;
m1 = bwdist(src2_mask);
m2 = bwdist(1 - hole_mask);
blend_alpha = m1 ./ (m1 + m2);

hole_mask(hole_mask ~= 0) = 1;
src2_mask(src2_mask ~= 0) = 1;



se = strel('disk',60);        
small_mask = imerode(src2_mask,se);
mapInd = find(repmat(small_mask,[1 1 3]));
im1_src(mapInd) = im2_src(mapInd);
im1 = im2_src;
mapInd = find(repmat(src2_mask,[1 1 3]));
im1(mapInd) = im2_src(mapInd);
mapInd = find(repmat(synth_mask,[1 1 3]));
im1(mapInd) = 0;

hole_mask = double(hole_mask);
src2_mask =  1 - double(src2_mask);
