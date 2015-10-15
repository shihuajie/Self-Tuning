function [output_im, output_NNF] = do_EM_iter_twoSRC(trg_im, src_im_ScalePyr, mask_Pyr, trg_mask, src_scales, options, num_em, scale_ind, big_size, org_im, org_mask,org_src1, org_src2, num_pyrs, alpha_mask,input_NNF_Pyrs)
global resize_method;
global save_intermediate;
global do_median;
global is_LAB;
global grad_weight;
global max_val_B;

num_src = 2;
cur_src = src_im_ScalePyr{scale_ind};
cur_alpha_mask = imresize(alpha_mask, [size(cur_src,1) size(cur_src,2)],resize_method);
num_channels = size(src_im_ScalePyr{scale_ind}, 3);
options.window_size = 0;
%src_mixed = [org_src1 org_src2];
cur_mask = mask_Pyr{scale_ind};
[cur_src, cur_mask_fixed] = fix_downsampled_mask(cur_src, cur_mask);
cur_org = cur_src;
cur_org(find(repmat(cur_mask_fixed, [1 1 num_channels]))) = -1000;

src1_im_PYR = make_pyramid_noscale(org_src1,zeros(size(org_src1,1),size(org_src1,2)), size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);
src2_im_PYR = make_pyramid_noscale(org_src2,zeros(size(org_src2,1),size(org_src2,2)), size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);
%src_mixed_PYR = make_pyramid_noscale(src_mixed,zeros(size(src_mixed,1),size(src_mixed,2)), size(trg_im,1)/big_size(1), options.max_scale, num_pyrs);%single(src_im_big);
src_im_PYR = {src1_im_PYR,src2_im_PYR};%,src_mixed_PYR};
alpha = cur_mask;
mask_cells = {repmat(cur_alpha_mask,[1 1 num_channels]), 1 - repmat(cur_alpha_mask,[1 1 num_channels])};
alpha(alpha < 0.05) = 0;
alpha_fixed = alpha;
alpha_fixed(alpha < 0.01) = 0;
%alpha_fixed(alpha > 0.99) = 1;
alpha_fixed_blend = repmat(cur_mask_fixed, [1 1 num_channels]);
%    alpha = max(alpha, double(cur_mask));
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


imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('C:\\Results\\Iterations\\Scale%03d_Iter%03d.png',(scale_ind), 0));

%options.EM_alpha = single(alpha_fixed(:,:,1));


for e = 1 : num_em
    trg_im_PYR.ims{1} = single(trg_im);
    trg_im_PYR.ims{1}(:,:,4) = 0;
    trg_im_PYR.ims{1}(:,:,5) = 0;
    trg_im_PYR.scales{1} = 1;

    for k = 1 : 2
        if ~isempty(prev_NNF_Pyrs{k})
            options.prev_map = prev_NNF_Pyrs{k};
        end
%         if(scale_ind > 1)
%             options.search_size = single(0.1*(1.1 - mask_cells{k}));
%         end
        
        opt = options;
        if e == 1
            opt.parallel_iter = 10;
            options.do_MinNNF = 4;
        end
        [cur_NNF GB] = snnMex(trg_im_PYR, src_im_PYR{k}, opt);
        GBArr{k} = GB;
        %cur_vote = svoteMedMex(trg_im_PYR, src_im_PYR{k}, cur_NNF, options, GB);
        cur_vote = svoteMex(trg_im_PYR, src_im_PYR{k}, cur_NNF, options, GB);
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
  
   %trg_im = alpha_fixed_blend .* trg_im + (1 - alpha_fixed_blend) .* cur_org;
    if(grad_weight)
        [dx dy]  = sGradMex(single(cur_org(:,:,1)));
%         if scale_ind >= 9
%             trg_im(:,:,4) = cur_mask_fixed .* trg_im(:,:,4) + (1 - cur_mask_fixed) .* dx;
%             trg_im(:,:,5) = cur_mask_fixed .* trg_im(:,:,5) + (1 - cur_mask_fixed) .* dy;
%         end
        %if(k == num_em)
        res = do_poisson_combination(trg_im(:,:,1), ones(size(trg_im(:,:,1))), grad_weight, trg_im(:,:,4), trg_im(:,:,5));
        %res = do_poisson_combination(trg_im(:,:,1), alpha_fixed, grad_weight, trg_im(:,:,4), trg_im(:,:,5));
        
        % imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('C:\\Results\\Iterations\\Scale%03d_Iter%03d_original.png',(scale_ind), k));
        trg_im(:,:,1) = res;%alpha_fixed(:,:,1) .* res + (1 - alpha_fixed(:,:,1)) .* cur_src(:,:,1);
        if scale_ind >= 9
            res(res > max_val_B) = max_val_B;
            trg_im(:,:,1) = res;
            %trg_im(:,:,1) = cur_mask_fixed .* trg_im(:,:,1) + (1 - cur_mask_fixed) .* cur_org(:,:,1);
            %trg_im(:,:,2) = cur_mask_fixed .* trg_im(:,:,2) + (1 - cur_mask_fixed) .* cur_org(:,:,2);
            %trg_im(:,:,3) = cur_mask_fixed .* trg_im(:,:,3) + (1 - cur_mask_fixed) .* cur_org(:,:,3);
        else
            res(res > 0.97 * max_val_B) = max_val_B;
            trg_im(:,:,1) = res;
        end
        %end
    end
    
     %trg_im = alpha_fixed_blend .* trg_im + (1 - alpha_fixed_blend) .* cur_src;
    
%     if e == 1
%        opt = options; 
%         opt.parallel_iter = 20;
%         opt.do_MinNNF = 0;
%         opt.sequential_iter = 0;
%         opt.patch_size = 8;
%         opt.min_scale = 0.99;%0.8;
%         opt.max_scale = 1.01;%1.5;
%         opt.min_angle = -pi;
%         opt.max_angle = pi;
%         
%         opt.min_Bias = single([-30    ; -30    ;   -30;      0;0]);
%         opt.max_Bias = single([30    ; 30    ;   30;      0;0]);
%         opt.min_Gain = single([0.9  ; 1.0     ;   1.0;      1;1]);
%         opt.max_Gain = single([1.1  ; 1.0   ;   1.0;    1;1]);
%         
%         opt.patch_size = 6;
%         opt.parallel_iter = 8;
%         [cur_NNF GB] = snnMex(trg_im_PYR, src_im_PYR{3}, opt);
%         opt.coh_weigh = 0;
%         cur_vote = svoteMex(trg_im_PYR, src_im_PYR{3}, cur_NNF, opt, GB);
%         out_im = cur_vote(:,:,1:num_channels);
%         out_im_arr{3} = out_im;
%         prev_NNF_Pyrs{3} = cur_NNF;
%         trg_im = out_im;%alpha .* out_im + (1 - alpha) .* cur_src;
%     end
    
     imwrite(lab2rgb(trg_im(:,:,1:3)), sprintf('C:\\Results\\Iterations\\Scale%03d_Iter%03d.png',(scale_ind), e));
    

    %imwrite(lab2rgb(out_im_arr{1}(:,:,1:3)), sprintf('C:\\Results\\Iterations\\A_Scale%03d_Iter%03d.png',(scale_ind), e));
    %imwrite(lab2rgb(out_im_arr{2}(:,:,1:3)), sprintf('C:\\Results\\Iterations\\B_Scale%03d_Iter%03d.png',(scale_ind), e));
    
    
    


end

 trg_im = alpha_fixed_blend .* trg_im + (1 - alpha_fixed_blend) .* cur_src;
output_im = trg_im;
output_NNF = prev_NNF_Pyrs;