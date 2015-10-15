function [output_im, output_NNF] = do_EM_iter_texture(trg_im, src_im_ScalePyr, mask_Pyr, trg_mask, src_scales, options, num_em, scale_ind, big_size, org_im, org_mask,org_src1, org_src2, num_pyrs, alpha_mask,input_NNF_Pyrs)

global out_directory;
global grad_weight;
global max_val_B;

num_src = 2;
cur_src = src_im_ScalePyr{scale_ind};
cur_alpha_mask = imresize(alpha_mask, [size(cur_src,1) size(cur_src,2)],'lanczos3');
num_channels = size(trg_im, 3); % 5
options.window_size = 0;

cur_mask = mask_Pyr{scale_ind};
[cur_src, cur_mask_fixed] = fix_downsampled_mask(cur_src, cur_mask);
cur_org = cur_src;
cur_org(find(repmat(cur_mask_fixed, [1 1 num_channels]))) = -1000;

src1_im_PYR = make_pyramid_noscale(org_src1,zeros(size(org_src1,1),size(org_src1,2)), size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);
src2_im_PYR = make_pyramid_noscale(org_src2,zeros(size(org_src2,1),size(org_src2,2)), size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);

src_im_PYR = {src1_im_PYR,src2_im_PYR};
alpha = cur_mask;
mask_cells = {repmat(cur_alpha_mask,[1 1 num_channels]), 1 - repmat(cur_alpha_mask,[1 1 num_channels])};
alpha(alpha < 0.05) = 0;
alpha_fixed = alpha;
alpha_fixed(alpha < 0.01) = 0;

alpha = repmat(alpha_fixed, [1 1 num_channels]);
weight_mask = single(mask2weights_hole(~alpha_fixed,1.3));

options.trg_weight = single(weight_mask);
options.EM_alpha = single(alpha(:,:,1));

prev_NNF_Pyrs = cell(num_src,1);
if ~isempty(input_NNF_Pyrs)
    for k = 1 : num_src
        pNN = single(input_NNF_Pyrs{k});
        NNF_NN = single(upscale_NNmap(pNN,size(trg_im,1), size(trg_im,2),size(src_im_PYR{k}.ims{1},1), size(src_im_PYR{k}.ims{1},2), size(src_im_PYR{k}.ims{1},1), size(src_im_PYR{k}.ims{1},2), options));
        prev_NNF_Pyrs{k} = NNF_NN;
    end
end


imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('%s\\Im%03d_Iter%03d.png', out_directory,(scale_ind), 0));
for e = 1 : num_em
    trg_im_PYR.ims{1} = single(trg_im);
    trg_im_PYR.ims{1}(:,:,4) = 0;
    trg_im_PYR.ims{1}(:,:,5) = 0;
    trg_im_PYR.scales{1} = 1;
    
    for k = 1 : 2
        if ~isempty(prev_NNF_Pyrs{k})
            options.prev_map = prev_NNF_Pyrs{k};
        end
        
        opt = options;
  
        [cur_NNF GB] = sNNMex(trg_im_PYR, src_im_PYR{k}, opt);

        % GBArr{k} = GB;
        cur_vote = svoteMex(trg_im_PYR, src_im_PYR{k}, cur_NNF, options, GB);
        out_im = cur_vote(:,:,1:num_channels);

        out_im_arr{k} = out_im;
        prev_NNF_Pyrs{k} = cur_NNF;
    end
    isMax = 1;%mod(e,2);
    
    trg_im = combine_two_images(out_im_arr{1}, out_im_arr{2},mask_cells{1},mask_cells{2},isMax);
    
    %trg_im = alpha_fixed_blend .* trg_im + (1 - alpha_fixed_blend) .* cur_org;
    if(grad_weight)
        [dx dy]  = sGradMex(single(cur_org(:,:,1)));
        %if scale_ind < 3
        trg_im(:,:,4) = cur_mask_fixed .* trg_im(:,:,4) + (1 - cur_mask_fixed) .* dx;
        trg_im(:,:,5) = cur_mask_fixed .* trg_im(:,:,5) + (1 - cur_mask_fixed) .* dy;
        %end
        %if(k == num_em)
        res = do_poisson_combination(trg_im(:,:,1), ones(size(trg_im(:,:,1))), grad_weight, trg_im(:,:,4), trg_im(:,:,5));
        trg_im(:,:,1) = res;%alpha_fixed(:,:,1) .* res + (1 - alpha_fixed(:,:,1)) .* cur_src(:,:,1);
        if scale_ind >= 9
            res(res > max_val_B) = max_val_B;
            trg_im(:,:,1) = res;
        else
            res(res > 0.95 * max_val_B) = max_val_B;
            trg_im(:,:,1) = res;
        end
        %end
    end
    
    
    imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('%s\\Im%03d_Iter%03d.png',out_directory,(scale_ind), e));
    
    
    
end

output_im = trg_im;
output_NNF = prev_NNF_Pyrs;