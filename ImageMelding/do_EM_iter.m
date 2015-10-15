function [output_im, output_NNF] = do_EM_iter(trg_im, src_im_ScalePyr, mask_Pyr, trg_mask, src_scales, options, num_em, cur_scale, big_size, org_im, org_mask, num_pyrs)
global grad_weight;
global out_directory;




src_im = src_im_ScalePyr{cur_scale};
src_im = cat(3,src_im, zeros(size(src_im,1),size(src_im,2),2));
num_channels = 5;

cur_src = src_im;
cur_mask = mask_Pyr{cur_scale};

[cur_src, cur_mask_fixed] = fix_downsampled_mask(cur_src, cur_mask);

if(grad_weight)%evaluate the gradients of boundaries

    [tmp1, tmp2] = sGradMex(single(cur_src(:,:, 1)));
    tmp1_f = imresize((tmp1),[size(org_im,1) size(org_im,2)], 'nearest');
    tmp2_f = imresize((tmp2),[size(org_im,1) size(org_im,2)], 'nearest');
    tmp1_f = imresize((tmp1_f),[size(cur_mask_fixed,1) size(cur_mask_fixed,2)], 'nearest');
    tmp2_f = imresize((tmp2_f),[size(cur_mask_fixed,1) size(cur_mask_fixed,2)], 'nearest');

    inds = find((cur_mask  > 0 & cur_mask < 1));
    tmp1(inds) = tmp1_f(inds);
    tmp2(inds) = tmp2_f(inds);

    cur_src(:,:,4) = tmp1_f;
    cur_src(:,:,5) = tmp2_f;

end

trg_im = imresize(trg_im, [size(cur_src,1) size(cur_src,2)] , 'lanczos3');

cur_trg_szX = size(trg_im,1);
cur_trg_szY = size(trg_im,2);

src_im_PYR = make_pyramid_noscale(org_im,org_mask,size(cur_src,1)/big_size(1), options.max_scale, num_pyrs);

alpha = cur_mask;

alpha(alpha < 0.01) = 0;

alpha = repmat(alpha, [1 1 num_channels]);
weight_mask = single(mask2weights_hole(~cur_mask_fixed,1.3));%gives more weights to regions closer to boundaries. For more detail see Wexler 2004

options.trg_weight = single(weight_mask);
options.EM_alpha = single(alpha(:,:,1));


if isfield(options, 'prev_map')
    pNN = single(options.prev_map);
    NNF_NN = single(upscale_NNmap(pNN,cur_trg_szX, cur_trg_szY,size(src_im_PYR.ims{1},1), size(src_im_PYR.ims{1},2), size(src_im_PYR.ims{1},1), size(src_im_PYR.ims{1},2), options));
    options = rmfield(options, 'prev_map');
    trg_im_PYR.ims{1} = single(trg_im);
    trg_im_PYR.scales{1} = 1;
    
    options.prev_map = NNF_NN;
end

imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('%s\\Im%03d_Iter%03d.png',out_directory,(cur_scale), 0));
alpha_fixed = alpha;
alpha_fixed(alpha < 0.01) = 0;
alpha_fixed(alpha > 0.99) = 1;

current_iter_num = options.parallel_iter;


for k = 1 : num_em
    opt = options;
    trg_im_PYR.ims{1} = single(trg_im);
    trg_im_PYR.scales{1} = 1;

    
    [cur_NNF, GB] = sNNMex(trg_im_PYR, src_im_PYR, opt);
 
    cur_vote = svoteMex(trg_im_PYR, src_im_PYR, cur_NNF, opt, GB);
    
    trg_im = cur_vote(:,:,1:num_channels);
    
    
    trg_im(:,:,1:3) = alpha_fixed(:,:,1:3) .* trg_im(:,:,1:3) + (1 - alpha_fixed(:,:,1:3)) .* cur_src(:,:,1:3);%put back the hole constain
    
    if(grad_weight)
        [dx, dy]  = sGradMex(trg_im(:,:,1));
        
        trg_im(:,:,4) = alpha_fixed(:,:,4) .* trg_im(:,:,4) + (1 - alpha_fixed(:,:,4)) .* dx;
        trg_im(:,:,5) = alpha_fixed(:,:,5) .* trg_im(:,:,5) + (1 - alpha_fixed(:,:,5)) .* dy;
        
        res = do_poisson_combination(trg_im(:,:,1), alpha_fixed(:,:,1), options.grad_weight(:,:,1), trg_im(:,:,4), trg_im(:,:,5));%solves screened Poisson equation
        
        
        trg_im(:,:,1) = res;
    end
    
    imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('%s\\Im%03d_Iter%03d.png',out_directory,(cur_scale), k));
    
    
    options.prev_map = cur_NNF;
end


output_im = alpha .* trg_im + (1 - alpha) .* cur_src;

output_NNF = cur_NNF;