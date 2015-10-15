function [trg_im, out_Pyr, mask_Pyr, out_scales, true_src_sizes] = generate_scale_pyramid_forEM( ...
    input_im, input_mask, start_scale, num_scales, is_logscale)

cur_scale = start_scale;
scale_step = (1 - start_scale) / (num_scales - 1);
im_szX = size(input_im,1); im_szY = size(input_im,2);

out_Pyr = cell(num_scales, 1);
mask_Pyr = cell(num_scales, 1);

true_src_sizes = cell(num_scales,1);
out_scales = zeros(num_scales,1);

logMin = log2(start_scale);

if(num_scales == 1)
    diff = 0;
    log_SZ = [0];
else
    diff = logMin / (num_scales - 1);
    log_SZ = 0 : diff : logMin;
end

sZs = flipdim(2 .^ log_SZ,2);
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
    
    src_im = imresize(input_im, [cur_szX cur_szY],'lanczos3');
    src_mask = imresize(input_mask,[cur_szX cur_szY],'lanczos3');
    if(k == 1)
        trg_im = src_im;
    end
    out_Pyr{k} = src_im;
    mask_Pyr{k} = src_mask;

end

