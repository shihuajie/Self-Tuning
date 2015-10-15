function out_pyr = make_pyramid_noscale(input_im, input_mask, min_scale, max_scale, num_scale)

global grad_weight;

logMin = log2(min_scale);
max_scale_global = min_scale + min_scale * max_scale;
logMax = log2(max_scale_global);
num_channels = size(input_im , 3);
if(num_scale == 1)
    diff = 0;
    log_SZ = [logMin];
else
    diff = (logMax - logMin) / (num_scale - 1);
    log_SZ = logMin : diff : logMax;
end

sZs = 2 .^ log_SZ;
ox  = size(input_im,1) ;
oy = size(input_im,2);
out_pyr_ims = cell(num_scale,1);
out_pyr.scales = cell(num_scale,1);
for k = 1 : num_scale 
    sZ = sZs(k);
    
    szn = round(sZ *[ox oy]);
    out_pyr.scales{k} = ox / szn(1);
    sZIM_mask = imresize(input_mask,szn,'lanczos3');
    
    sZIM = imresize(input_im,szn,'lanczos3');
    sZIM_mask_b = imresize(sZIM_mask, [ox oy], 'lanczos3');
    sZIM_b = imresize(sZIM, [ox oy], 'lanczos3');
    
    [sZIM_b, sZIM_mask_b] = fix_downsampled_mask(sZIM_b, sZIM_mask_b);
    sZIM_b(find(repmat(sZIM_mask_b,[1 1 num_channels]))) = -1000;
  
    out_pyr_ims{k} = single(sZIM_b);
end

out_pyr.ims = out_pyr_ims;
