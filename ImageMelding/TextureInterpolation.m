% Image Melding 
% MATLAB Scripts and Mex Executables Version 1.0 (2012-03-12)
% 
%   By Soheil Darabi
%   Copyright Adobe Systems Inc
%   Licensed by Adobe for noncommercial research use only.

clc;
clear options;
close all;
path(path, '.\Mexfiles');

global grad_weight;%lambda value in the paper
global out_directory;%directory name for saving intermediate files

input_file_name = 'Bricks_Plastic.png';
out_directory = '.\Results';
grad_weight = 2.0;

 im1_src = im2double(imread('Plastic.png')); 
 im2_src = im2double(imread( 'Bricks.png'));
[im1 im1_mask invalid_mask] = read_mask_image(input_file_name,1,0.8); %reads the input image and wherever the color is magenta it assings it to be a hole!

 alpha_mask = generate_linear_blend(im1_mask);    

options.parallel_iter =6; %number of random search in Generalized PatchMatch (GPM)
options.patch_size = 10;%patchsize
options.min_scale = 0.9;%minimum scale in GPM 
options.max_scale = 1.3;%maximum scale in GPM 
options.min_relscale = 0.9;%minimum relative scale in GPM 
options.max_relscale = 1.1;%maximum relative scale in GPM 
options.min_angle = -pi/4;%minimum rotation
options.max_angle = pi/4;%maximum rotation
options.search_reflection = 1;%indicates if the algorithm should search for reflection.  If you turn it on you need more parallel_iter!
options.window_size = 0;%Specifies the local search
options.random_init = 1;%if it is turned to zero, GPM initializes the locations to be identity
options.do_random_search = 1;%if it is turned to zero, algorithm just does random search

options.min_Bias = single([-10   ; -50   ;   -50;      0; 0]);%minimum bias in Lab space.  The last two channels are related to gradients.  As gradients can be drived from colors, if you leave it as it's the algorithm automatically uses gain of color channels.
options.max_Bias = single([10 ; 50 ;  50; 0;0]);%maximum bias 
options.min_Gain = single([0.9  ; 0.9     ;   0.9;      1;1]);%minimum gain
options.max_Gain = single([1.5  ; 1.5   ;   1.5;    1;1]);%maximum gain
init_im = random_merge(im1_src,im2_src,alpha_mask,options.patch_size);%the initial guess is created by selecting patches from each source randomly.  The probablity of the selcetion is based on how close the synthesis region is to each of the sources.

options.coh_weigh = 1;%if set to 1 the algorithm gives more weight to the coherent regions.
options.normalize_weight = 0;%if set to 0 the algorithm ignores gain and bias adjustments
options.num_threads = 12;%number of threads
options.do_MinNNF = 2;%the algorithm does this number of random searches after merges previous NNF and current one.  Please refer to original PM for more details
options.grad_weight = grad_weight;
num_pyrs = 10;%number of scaled pyramids algorithm needs for synthesis

big_size =  [size(im1,1), size(im1,2)];

start_scale = 35 / min(size(im1,1), size(im1,2));%the coasest scale is 35

num_scales = 10;%number of scales for EM iterations
num_em = 30;%initial number of EM iterations at coarsest scale
decrease_factor = 3;%when we go to finer resolutions we reduce this number from number of EM iterations
min_scale = 5;%the minimum number of EM iterations 

C = makecform('srgb2lab');
im1 = applycform(im1, C);%takes the image to Lab space
im1 = cat(3,im1,zeros(size(im1,1),size(im1,2),2));%add two channels for gradients.
im1_src = applycform(im1_src, C);%takes the image to Lab space
im1_src = cat(3,im1_src,zeros(size(im1_src,1),size(im1_src,2),2));%add two channels for gradients.
im2_src = applycform(im2_src, C);%takes the image to Lab space
im2_src = cat(3,im2_src,zeros(size(im2_src,1),size(im2_src,2),2));%add two channels for gradients.


%options.vote_weight = single(fspecial('gaussian', options.patch_size,options.patch_size/1));%single(ones(options.patch_size,options.patch_size));
options.vote_weight = single(ones(options.patch_size,options.patch_size));%instead of box filter in vote step you can apply other filters such as Gaussian

im_szX = size(im1,1); im_szY = size(im1,2);
is_log_upscale = 1;

[trg_im,  src_Pyr, mask_Pyr, src_scales, true_src_sizes] = generate_scale_pyramid_forEM(im1, im1_mask, start_scale, num_scales,is_log_upscale);%generates all the source pyramids for every EM scale

trg_mask = mask_Pyr{1};
[trg_im,trg_mask] = fix_downsampled_mask(trg_im, trg_mask); %when downsample an image with a hole, the values of the boundaries get affected by both non-valid regions and valid regions.  This function correct those values.

trg_im(:,:,1:3) = imresize(applycform(init_im, C),[size(trg_im,1),size(trg_im,2)]);%applycform(tmp_im, C);

global max_val_B;
max_val_B =  max(max(trg_im(:,:,1)));%to control the colors don't get to brigh in synthesis I clamp values bigger than propertion of the maximum value.
cur_NN = [];;

for k = 1 : num_scales%loop over scales
    
    if k ~= 1
        
        trg_mask = imresize(trg_mask, true_src_sizes{k},'lanczos3');
        trg_im = imresize(cur_im, true_src_sizes{k},'lanczos3');
    end
    
    opt = options;
    if k > 3%for acceleration I avoid doing random search after a specific scale
        opt.do_random_search = 0;
        opt.do_MinNNF = 0;
    end
    size(im1, 1)
    size(im1, 2)
    big_size =  [size(im1,1), size(im1,2)]
    [cur_im cur_NN] = do_EM_iter_texture(trg_im, src_Pyr, mask_Pyr, trg_mask, src_scales, opt, num_em, k,big_size, im1, im1_mask, im1_src, im2_src, num_pyrs,alpha_mask,cur_NN);
    
    num_em = max(min_scale, num_em - decrease_factor);
end
