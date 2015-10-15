function [output_im, output_NNF] = do_EM_iter_twoSRC_cloning(trg_im, src_im_ScalePyr, mask_Pyr, trg_mask, src_scales, options, num_em, scale_ind, big_size, org_im, org_mask,org_src1, org_src2, num_pyrs, alpha_mask,input_NNF_Pyrs, hole_mask, src2_mask)

global grad_weight;
global max_val_B;
global out_directory;

num_src = 2;
non_valid_width =  (1 - 0.0) / 2;
cur_src = src_im_ScalePyr{scale_ind};
cur_alpha_mask = imresize(alpha_mask, [size(cur_src,1) size(cur_src,2)],'lanczos3');
num_channels = size(src_im_ScalePyr{scale_ind}, 3);


cur_mask = mask_Pyr{scale_ind};
[cur_src, cur_mask_fixed] = fix_downsampled_mask(cur_src, cur_mask);

cur_org = cur_src;
cur_org(find(repmat(cur_mask_fixed, [1 1 num_channels]))) = -1000;
cur_size = [size(trg_im,1) size(trg_im,2)];

non_valid_mask = zeros(size(alpha_mask));
non_valid_mask(alpha_mask > non_valid_width & alpha_mask < 1 - non_valid_width) = 1;

src1_im_PYR = make_pyramid_noscale_known(org_src1,non_valid_mask, size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);
src2_im_PYR = make_pyramid_noscale_known(org_src2,non_valid_mask, size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);
%src_mixed_PYR = make_pyramid_noscale(src_mixed,zeros(size(src_mixed,1),size(src_mixed,2)), size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);
src_im_PYR = {src1_im_PYR,src2_im_PYR};%,src_mixed_PYR};
alpha = cur_mask;
cur_alpha_mask = sqrt(abs(cur_alpha_mask));
mask_cells = {repmat((cur_alpha_mask),[1 1 num_channels]), 1 - repmat((cur_alpha_mask),[1 1 num_channels])};
alpha(alpha < 0.05) = 0;
alpha_fixed = alpha;
alpha_fixed(alpha < 0.01) = 0;

im1_mask_blend = repmat(imresize(hole_mask,cur_size,'lanczos3'),[1 1 num_channels]);
im2_mask_blend = repmat(imresize(src2_mask,cur_size,'lanczos3'),[1 1 num_channels]);

im1_mask_blend_fixed = im1_mask_blend;
im1_mask_blend_fixed(im1_mask_blend_fixed < 0.5) = 0;
im1_mask_blend_fixed(im1_mask_blend_fixed ~= 0) = 1;

im2_mask_blend_fixed = im2_mask_blend;
im2_mask_blend_fixed(im2_mask_blend_fixed < 0.5) = 0;
im2_mask_blend_fixed(im2_mask_blend_fixed ~= 0) = 1;

im1_blend_src = imresize(org_src1,cur_size,'lanczos3');
im2_blend_src = imresize(org_src2,cur_size,'lanczos3');
[im1_blend_src(:,:,4) im1_blend_src(:,:,5)]  = sGradMex(single(im1_blend_src(:,:,1)));
[im2_blend_src(:,:,4) im2_blend_src(:,:,5)]  = sGradMex(single(im2_blend_src(:,:,1)));
grad1_mask = imfilter(im1_mask_blend,[0.25 .25;0.25 .25]);
grad2_mask = imfilter(im2_mask_blend,[0.25 .25;0.25 .25]);
%grad1_mask(grad1_mask ~= 0) =  1;
%grad2_mask(grad2_mask ~= 0) =  1;
%alpha_fixed(alpha > 0.99) = 1;
alpha_fixed_blend = repmat(cur_mask_fixed, [1 1 num_channels]);
%    alpha = max(alpha, double(cur_mask));
alpha = repmat(alpha_fixed, [1 1 num_channels]);
alpha_weight = 1.3;
%weight_mask = single(mask2weights_hole(~alpha_fixed,1.01));
weight_mask{1} = single(mask2weights_hole(~im1_mask_blend_fixed,1.3));
weight_mask{2} = single(mask2weights_hole(~im2_mask_blend_fixed,1.3));

options.EM_alpha = single(alpha(:,:,1));

prev_NNF_Pyrs = cell(num_src,1);
if ~isempty(input_NNF_Pyrs)
    for k = 1 : num_src
        pNN = single(input_NNF_Pyrs{k});
        NNF_NN = single(upscale_NNmap(pNN,size(trg_im,1), size(trg_im,2),size(src_im_PYR{k}.ims{1},1), size(src_im_PYR{k}.ims{1},2), size(src_im_PYR{k}.ims{1},1), size(src_im_PYR{k}.ims{1},2), options));
        prev_NNF_Pyrs{k} = NNF_NN;
    end
end


%options.EM_alpha = single(alpha_fixed(:,:,1));
trg_im(:,:,1:3) = im1_mask_blend(:,:,1:3) .* trg_im(:,:,1:3) + (1 - im1_mask_blend(:,:,1:3)) .* im1_blend_src(:,:,1:3);
%trg_im(:,:,1:3) = im2_mask_blend(:,:,1:3) .* trg_im(:,:,1:3) + (1 - im2_mask_blend(:,:,1:3)) .* im2_blend_src(:,:,1:3);

trg_im(:,:,4:5) = grad1_mask(:,:,4:5) .* trg_im(:,:,4:5) + (1 - grad1_mask(:,:,4:5)) .* im1_blend_src(:,:,4:5);
trg_im(:,:,4:5) = grad2_mask(:,:,4:5) .* trg_im(:,:,4:5) + (1 - grad2_mask(:,:,4:5)) .* im2_blend_src(:,:,4:5);


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
        options.trg_weight = single(weight_mask{k}(:,:,1));
        %         if(scale_ind > 1)
        %             options.search_size = single(0.1*(1.1 - mask_cells{k}));
        %         end
        
        opt = options;
        if e == 1
            if scale_ind < 4
                opt.do_MinNNF = 4;
                opt.parallel_iter = 10;
                
                
                %opt.search_size = single(0.1*ones(cur_size));
                
            else
                %             opt.parallel_iter = 20;
                %             opt.do_MinNNF = 4;
                opt.do_MinNNF = 0;
                opt.do_random_search = 0;
                
                
                %opt.search_size = single(0.1*ones(cur_size));
            end
        end
        [cur_NNF GB] = sNNMex(trg_im_PYR, src_im_PYR{k}, opt);
        GBArr{k} = GB;
        %cur_vote = svoteMedMex(trg_im_PYR, src_im_PYR{k}, cur_NNF, options, GB);
        cur_vote = svoteMex(trg_im_PYR, src_im_PYR{k}, cur_NNF, opt, GB);
        out_im = cur_vote(:,:,1:num_channels);
        %[out_im(:,:,4) out_im(:,:,5)]  = sgradMex(out_im(:,:,1));
        
        %out_im = ( mask_cells{k}) .* out_im + ( 1 - mask_cells{k}) .* cur_src;
        
        %res = do_poisson_combination(out_im(:,:,1), double(hole_alpha(:,:,1)), double(options.grad_weight), double(out_im(:,:,4)), double(out_im(:,:,5)));
        %out_im(:,:,1) = res;
        out_im_arr{k} = out_im;
        prev_NNF_Pyrs{k} = cur_NNF;
    end
    
    %trg_im = ( mask_cells{1} .* out_im_arr{1} +  mask_cells{2} .* out_im_arr{2}) ./ (mask_cells{1}+mask_cells{2});
    %trg_im(:,:,4:5) =  max(out_im_arr{1}(:,:,4:5),out_im_arr{2}(:,:,4:5));
    isMax = 1;%mod(e,2);
    
    trg_im = combine_two_images(out_im_arr{1}, out_im_arr{2},mask_cells{1},mask_cells{2},isMax);
    
    trg_im(:,:,1:3) = im1_mask_blend(:,:,1:3) .* trg_im(:,:,1:3) + (1 - im1_mask_blend(:,:,1:3)) .* im1_blend_src(:,:,1:3);
    %trg_im(:,:,1:3) = im2_mask_blend(:,:,1:3) .* trg_im(:,:,1:3) + (1 - im2_mask_blend(:,:,1:3)) .* im2_blend_src(:,:,1:3);
    
    trg_im(:,:,4:5) = grad1_mask(:,:,4:5) .* trg_im(:,:,4:5) + (1 - grad1_mask(:,:,4:5)) .* im1_blend_src(:,:,4:5);
    trg_im(:,:,4:5) = grad2_mask(:,:,4:5) .* trg_im(:,:,4:5) + (1 - grad2_mask(:,:,4:5)) .* im2_blend_src(:,:,4:5);
    
    
    res = do_poisson_combination(trg_im(:,:,1), alpha_fixed, grad_weight, trg_im(:,:,4), trg_im(:,:,5));
    res(res > max_val_B) = max_val_B;
    % imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('C:\\Results\\Iterations\\Scale%03d_Iter%03d_original.png',(scale_ind), k));
    trg_im(:,:,1) = res;%alpha_fixed(:,:,1) .* res + (1 - alpha_fixed(:,:,1)) .* cur_src(:,:,1);
    trg_im(:,:,1:3) = im1_mask_blend(:,:,1:3) .* trg_im(:,:,1:3) + (1 - im1_mask_blend(:,:,1:3)) .* im1_blend_src(:,:,1:3);
    
    imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('%s\\Im%03d_Iter%03d.png',out_directory,(scale_ind), e));
end







%trg_im = alpha_fixed_blend .* trg_im + (1 - alpha_fixed_blend) .* cur_src;
output_im = trg_im;
output_NNF = prev_NNF_Pyrs;