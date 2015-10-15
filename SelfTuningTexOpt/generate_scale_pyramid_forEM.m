function [trg_im, out_Pyr, mask_Pyr, out_scales, true_src_sizes] = generate_scale_pyramid_forEM( ...
    input_im, input_mask, start_scale, num_scales, is_logscale, ...
    multires, patch_size, multires_weight)
  
if nargin < 6
   multires = 0;
elseif nargin < 7
   error('Multi-resolution requires a patch size!');
end
if multires && nargin < 8
   multires_weight = 1; 
end

% scale numbers
cur_scale = start_scale;
scale_step = (1 - start_scale) / (num_scales - 1);

% size of the input
im_szX = size(input_im,1); 
im_szY = size(input_im,2);

% result cells (img+mask)
out_Pyr = cell(num_scales, 1 + multires);
mask_Pyr = cell(num_scales, 1);

% result sizes and scales
true_src_sizes = cell(num_scales,1);
out_scales = zeros(num_scales,1);

% pyramid scales
if(num_scales == 1)
    log_SZ = [0];
else
    logMin = log2(start_scale);
    diff = logMin / (num_scales - 1);
    log_SZ = 0 : diff : logMin;
end
sZs = flipdim(2 .^ log_SZ, 2);

% generate the pyramid
for k = 1 : num_scales
    if is_logscale
        cur_scale = sZs(k);
        cur_szX = round(cur_scale * im_szX);
        cur_szY = round(cur_scale * im_szY);
    else
        cur_szX = round(cur_scale * im_szX);
        cur_szY = round(cur_scale * im_szY);
        cur_scale = cur_scale + scale_step;
    end
    true_src_sizes{k} = [cur_szX cur_szY];
    out_scales(k) = cur_scale;
    
    % resized image
    %% modified by huajie 2015-9-15
    label_im = input_im(:,:,8);
    src_im = imresize(input_im, [cur_szX cur_szY], 'lanczos3');
    label_im = imresize(label_im, [cur_szX cur_szY], 'nearest');
    src_im(:, :, 8) = label_im;
    %% end
%     src_im = imresize(input_im, [cur_szX cur_szY], 'lanczos3');
    out_Pyr{k, 1} = src_im;
    
    % resized mask
    if ~isempty(input_mask)
        src_mask = imresize(input_mask, [cur_szX cur_szY], 'bilinear');
    else
        src_mask = [];
    end
    mask_Pyr{k, 1} = src_mask;
    
    % the starting image
    if(k == 1)
        trg_im = src_im;
    end
    
    % generate the texture layers for multi-resolution
    if multires >= 1
        out_Pyr(k, :) = create_multires(src_im, multires, patch_size, multires_weight);
    end
end
