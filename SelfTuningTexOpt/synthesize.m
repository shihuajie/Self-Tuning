function [ res_im, res_nnf ] = synthesize( src_im, org_im, org_mask, options )
%SYNTHESIZE - synthesize a full texture
%
% INPUT
%   - src_im    the exemplar image (L, a, b, Gx, Gy, F ...)
%   - org_im    the original result (L, a, b, Gx, Gy, F ...)
%   - org_mask  the original mask (single)
%   - options   the diverse synthesis options
%
% OUTPUT
%   - res_im    the synthesized image
%   - res_nnf   the last resulting nnf
%
    global out_directory;
    
    %% Loading mat environment
%     % load environment
%     env_file = trg_im;
%     load(env_file);
%     % restore rng state
%     rng(rng_state);
%     % move to safe directory
%     if ~exist(out_directory, 'dir')
%         out_directory = fileparts(env_file);
%     end
%     num = 0;
%     while exist([out_directory filesep 'tmpEnv' num2str(num)], 'file')
%         num = num + 1;
%     end
%     out_directory = [out_directory filesep 'tmpEnv' num2str(num)];
%     mkdir(out_directory); % create the environment directory
% 
%     % info
%     fprintf('... sandboxed in %s.\n... resuming at scale %d, iter %d.\n\n', out_directory, k, iter);
        
    %% Image pyramids
    tx_szX = size(src_im, 1);
    tx_szY = size(src_im, 2);
    % fixed scale
    if isfield(options, 'fixed_scale')
        scale_factor = options.fixed_scale;
        side = min([tx_szX, tx_szY]);
        num_scales = 0;
        while side >= options.coarsest_scale
            num_scales = num_scales + 1;
            side = side / scale_factor;
        end
        options.num_scales = num_scales;
        fprintf('Scale factor=%f, min=%f, num=%f\n', scale_factor, side, num_scales);
    end
    % the starting scale brings the pictures sides to at least options.coarsest_scale
    %% added by huajie 2015-9-17
    if(min([tx_szX, tx_szY]) > 1500) 
        options.coarsest_scale = 50;
    end
    %% end
    start_scale = options.coarsest_scale / min([tx_szX, tx_szY]);

    % original image pyramid
    [trg_im, org_Pyr, msk_Pyr, ~, true_trg_sizes] = generate_scale_pyramid_forEM( ...
      org_im, org_mask, start_scale, options.num_scales, options.is_log_upscale);
    multires = 0;
    multires_weight = 1;
    if isfield(options, 'multires')
        multires = options.multires;
        % if it's an array, it changes with each scale index
        if numel(multires) > 1
            multires = max(multires);
        else
            options.multires = [0, 0, ceil((3:options.num_scales) * multires / options.num_scales)];
        end
        % we need a weight!
        if isfield(options, 'multires_weight')
            multires_weight = options.multires_weight;
            if numel(multires_weight) == 1 && multires > 1
                options.multires_weight = ones(1, multires) * multires_weight;
            end
        else
            options.multires_weight = ones(1, multires);
        end
    else
        options.multires = 0;
    end

    % exemplar pyramid
    [~, tex_Pyr, ~, ~, true_tex_sizes] = generate_scale_pyramid_forEM( ...
      src_im, [], start_scale, options.num_scales, options.is_log_upscale, ...
      multires, options.patch_size, multires_weight);
    % correct the gradient
    for k = 1:numel(tex_Pyr)
        img = tex_Pyr{k, 1};
        [gx, gy] = gradient1(img(:, :, 1));
        img(:, :, options.grad_channels) = cat(3, gx, gy);
        tex_Pyr{k} = img;
    end
    
%    imwrite(lab2rgb(src_start_im(:,:,1:3)), [out_directory 'src_start_im.png']);
     
    
%% added by huajie 2015-9-9
%% start replace the guidance feature channel
%      % trg_gc_py
%      trg_gc_filepath = strcat(out_directory, '\trg_gc.png');
%      trg_gc_im = im2double(imread(trg_gc_filepath));
%      trg_gc_im = improc(trg_gc_im(:,:,1), 'otsu,rev,sbwdist,n');
%       trg_gc_filepath = strcat(out_directory, '\test_bin_bwdis.png');
%       imwrite(trg_gc_im, trg_gc_filepath);
% 
%      [trg_gc_feature, trg_gc_Pyr, ~, ~, true_trg_sizes] = generate_scale_pyramid_forEM( ...
%        trg_gc_im, org_mask, start_scale, options.num_scales, options.is_log_upscale);
%   
%        for k = 1:numel(trg_gc_Pyr)
%         img = trg_gc_Pyr{k, 1};
%         org_im = org_Pyr{k, 1};
%         org_im(:,:,6) = img;
%         org_Pyr{k,1} = org_im;
% %         trg_gc_image = mat2gray(img(:,:));
% %         imwrite(trg_gc_image, [out_directory filesep ...
% %             sprintf('trg_gc_s%02d.png', k) ]);
%        end
%% end
%% start build the target guidance channel pyramid
%% start replace the guidance feature channel
%      % added by huajie 2015-9-15
%      str = out_directory;
%      index_dir = findstr(str,'\');
%      str_temp = str(1:index_dir(end)-1);
%      src_gc_filepath = strcat(str_temp, '\src_gc.png');
%      src_gc_im = double(imread(src_gc_filepath));
%      trg_gc_filepath = strcat(str_temp, '\trg_gc.png');
%      trg_gc_im = double(imread(trg_gc_filepath));
% 
%       [src_gc_start, src_gc_Pyr, ~, ~, true_tex_sizes] = generate_scale_pyramid_forEM( ...
%       src_gc_im, [], start_scale, options.num_scales, options.is_log_upscale, ...
%       multires, options.patch_size, multires_weight);
%       for k = 1:numel(src_gc_Pyr)
%         img = src_gc_Pyr{k, 1};
%         tex_im = tex_Pyr{k, 1};
%         tex_im(:,:,8) = rgb2label(img)
%         tex_Pyr{k,1} = tex_im;
%         src_gc_image = mat2gray(tex_im(:,:,8));
%         imwrite(src_gc_image, [out_directory filesep ...
%              sprintf('src_gc_s%02d.png', k) ]);
%       end
%      
%      [trg_gc_start, trg_gc_Pyr, ~, ~, true_trg_sizes] = generate_scale_pyramid_forEM( ...
%        trg_gc_im, org_mask, start_scale, options.num_scales, options.is_log_upscale);
%      trg_im = cat(3, trg_im, rgb2label(trg_gc_start));
%      for k = 1:numel(trg_gc_Pyr)
%         img = trg_gc_Pyr{k, 1};
%         org_im = org_Pyr{k, 1};
%         y = rgb2label(img)
%         cat(3, org_im, y);
%         org_im(:,:,8) = y;
%         org_Pyr{k,1} = org_im;
%         trg_gc_image = mat2gray(org_im(:,:,8));
%         imwrite(trg_gc_image, [out_directory filesep ...
%              sprintf('trg_gc_s%02d.png', k) ]);
%      end
%% end
    %% Data parameters
    if ~isfield(options, 'max_brightness')
        % by default, we clamp brighter pixels to the maximum of the source
        options.max_brightness = max(max(src_im(:, :, 1)));
    end

    %% Scale loop and EM iterations
    cur_im = trg_im;
    convData = [];
    
    % use a previous nnf if there is one
    scramble_nnf = 0;
    if isfield(options, 'lastNNF')
        cur_nnf = load_mat(options.lastNNF);
        cur_nnf = rescale_NNF(cur_nnf, true_trg_sizes{1}, true_tex_sizes{1}, true_tex_sizes{end}, options);
        if ~isempty(org_mask)
            scramble_nnf = find(org_mask < 0, 1);
        end
    else
        cur_nnf = []; % no previous nnf to start with
    end
    rev_nnf = [];
    if strcmp(options.vote_method, 'bidir_sim')
        options.reverse_nnf = 1;
    end

    % vectorizable parameters depending on the scale level
    vectorizable = { 'iterations', 'multires', 'until_convergence', ...
      'coh_rand_search', 'max_rand_search', 'window_size', ...
      'comp_penalty', 'comp_weight', 'comp_exponent', ...
      'comp_abs_threshold', 'comp_rel_threshold'
    };

    % counters
    iter = 1; % iteration number
    k = 1; % scale number
    j = 0; % previous scale (before k)
    
    % until parameters to look at
    until_values = {};
    option_names = fieldnames(options);
    for f = 1:numel(option_names)
        opt_name = option_names{f};
        until_idx = strfind(opt_name, '_until');
        if ~isempty(until_idx) && until_idx == numel(opt_name) - 5 % 6 for _until, + 1
            target_name = opt_name(1:(until_idx - 1));
            if ~isfield(options, target_name)
                continue
            end
            until_values{end + 1} = target_name;
            last_scale = options.(opt_name);
            if isfield(options, 'fixed_scale')
                last_scale = round(last_scale * options.num_scales / 10);
            end
            until_values{end + 1} = last_scale;
            if isfield(options, [target_name '_then'])
                until_values{end + 1} = options.([target_name '_then']);
            else
                until_values{end + 1} = 0;
            end
        end
    end
    
    start_time = tic; % dependent of execution instance, not variable environment!
    start_em = options.num_em;
    % Loop
    while k <= options.num_scales

        % disable elements past their _until value
        for u = 1:3:numel(until_values)
            until_name = until_values{u};
            until_value = until_values{u + 1};
            until_then = until_values{u + 2};
            if k > until_value && isfield(options, until_name)
                if until_then
                    options.(until_name) = until_then;
                    fprintf('Replacing %s\n', until_name);
                elseif numel(options.(until_name)) > 1
                    options = rmfield(options, until_name);
                    fprintf('Removing %s option.\n', until_name);
                else
                    options.(until_name) = 0;
                    fprintf('Setting %s to zero.\n', until_name);
                end
            end
        end
      
        %% Apply the scale change if needed
        if j > 0
            fprintf('-----\nScaling from s%d (%d x %d) to s%d (%d x %d)\n-----\n', j, true_trg_sizes{j}(1), true_trg_sizes{j}(2), k, true_trg_sizes{k}(1), true_trg_sizes{k}(2));
            [cur_im, cur_nnf, ~] = rescale(cur_im, cur_nnf, [], ...
                true_trg_sizes, true_tex_sizes, k, j, options);
            if isfield(options, 'reverse_nnf') && numel(rev_nnf) > 1
                [~, rev_nnf, ~] = rescale(tex_Pyr{j, 1}, rev_nnf, [], ...
                    true_tex_sizes, true_trg_sizes, k, j, options);
            end
        end
        if isfield(options, 'reverse_nnf')
            cur_nnf = { cur_nnf, rev_nnf };
        end
        fprintf('\n\nScale #%02d (it=%d)\n\n', k, iter);
        
        %% Current mask
        cur_mask = msk_Pyr{k};
        % maybe save mask
        if options.saveMask
            mask = min(max(cur_mask > 0, 0), 1);
            if isempty(mask)
                P = options.patch_size;
                mask = zeros(size(cur_im, 1) - P + 1, size(cur_im, 2) - P + 1);
            end
            imwrite(mask, [out_directory filesep ...
              sprintf('mask-%03d_s%02d.png', iter, k) ]);   
        end

        %% The local options
        opt = options;
        
        % scaled information
        if k < options.num_scales
            % correct the generator and jitter size
            down_scale_factor = true_trg_sizes{k} ./ true_trg_sizes{options.num_scales};
            if isfield(options, 'generators')
                G = options.generators;
                G(1:2, 1:2) = G(1:2, 1:2) * down_scale_factor(1);
                opt.generators = G;
            end
            if isfield(options, 'align_jitter')
                opt.align_jitter = options.align_jitter * down_scale_factor(1);
            end
        end
        if k > 1
            up_scale_factor = true_trg_sizes{k} ./ true_trg_sizes{1};
            if isfield(options, 'segsize_threshold')
                opt.segsize_threshold = up_scale_factor(1) * options.segsize_threshold; 
            end
        end

        % vectors to numbers
        for v = 1:numel(vectorizable)
            vec_name = vectorizable{v};
            if isfield(opt, vec_name)
                value = opt.(vec_name);
                if numel(value) >= options.num_scales
                    opt.(vec_name) = value(k); % keep only the scale's value
                end
            end
        end
        
        % whether to scramble the first nnf
        if scramble_nnf
            % we should scramble the nnf
            opt.scramble_nnf = 1;
            scramble_nnf = 0; % not the next time
        end
        
        % scale and iteration numbers
        opt.last_scale = k == options.num_scales;
        num_em = max(options.min_em, start_em - (k - 1) * options.decrease_factor);
%% add by huajie 2015-9-9
        cur_im(:, :, 8) = org_Pyr{k,1}(:, :, 8);
 %% end 
        %% The EM iterations
        if options.saveConvData
            [cur_im, cur_nnf, ~, iter, cur_conv] = do_EM_iteration( ...
                cur_im, tex_Pyr, org_Pyr{k}, ...
                iter, num_em, k, cur_nnf, cur_mask, opt);
            convData = [convData, cur_conv];
        else
            [cur_im, cur_nnf, ~, iter] = do_EM_iteration(cur_im, tex_Pyr, org_Pyr{k}, ...
                iter, num_em, k, cur_nnf, cur_mask, opt);
        end
        
        % unpack nnf
        if iscell(cur_nnf)
            rev_nnf = cur_nnf{2};
            options.reverse_nnf = rev_nnf;
            cur_nnf = cur_nnf{1};
        end

        %% Scale change
        j = k;
        k = k + 1;

        % saving the last environment
        if options.saveEnv && k > options.num_scales
            % store the rng state
            rng_state = rng();
            % save environment
            save([out_directory filesep sprintf('env_s%02d_i%03d.mat', ...
                k, iter)]);
        end
    %% added by huajie 2015-8-31
    % save result for posterity
%     img = tex_Pyr{k-1, 1};
%     src_mid_im = lab2rgb(img(:,:,1:3));
%     imwrite(src_mid_im, [out_directory filesep ...
%              sprintf('src_mid_im%02d.png', k-1) ]);
    trg_mid_im = lab2rgb(cur_im(:, :, 1:3));
    imwrite(trg_mid_im, [out_directory filesep ...
             sprintf('trg_mid_im%02d.png', j) ]);

     trg_mid_feature_im = im2double(mat2gray(cur_im(:,:,8)));
     imwrite(trg_mid_feature_im, [out_directory filesep ...
              sprintf('cur_mid_label_im%02d.png', j) ]);
     % draw the rotation angle arrow
     nnf_angle_im = drawAngleArrow(cur_nnf(:,:,3));
     imwrite(nnf_angle_im, [out_directory filesep ...
              sprintf('nnf_angle_im%02d.png', j) ]);       
    %% end
    end

    %% Output
    res_im = cur_im;
    res_nnf = cur_nnf;

    %% Should we save the NNF matrix?

    if options.saveConvData
        save([out_directory filesep 'convData.mat'], 'convData');
    end
    
    if options.saveLastNNF
      save([out_directory filesep 'lastNNF.mat'], 'res_nnf');
    end
    
    % elapsed time
    elapsed_time = toc(start_time)
    save_data(elapsed_time, 'elapsed');
        
    
    % touch the "done" file
    fclose(fopen(fullfile(out_directory, 'done'), 'w'));

end

%% Rescaling (general)
function [new_im, new_nnf, new_mask] = rescale(im, nnf, mask, trg_sizes, tex_sizes, k, j, options)
    
    % resizing the image
    if k < j
        method = 'lanczos3'; % use lanczos(3) when downscaling to increase detail
    else
        method = 'bilinear'; % use linear interpolation when upsampling to be close to nnf
    end
    new_im = imresize(im, trg_sizes{k}, method);
    
    % rescaling the nnf
    trg_sz = trg_sizes{k};
    prev_tex_sz = tex_sizes{j};
    tex_sz = tex_sizes{k};
    new_nnf = rescale_NNF(nnf, trg_sz, tex_sz, prev_tex_sz, options);
    
    % create the mask
    if numel(mask) > 1
        new_mask = zeros(trg_sizes{j}, 'single');
        S = options.patch_size - 1;
        new_mask(1:(end - S), 1:(end - S)) = mask;
        % downscaling
        new_mask = imresize(new_mask, trg_sizes{k}, 'bilinear'); % TODO should we clean it?
        % cropping the mask
        new_mask = new_mask(1:(end - S), 1:(end - S));
        
        % cleaning negative if we are upsampling
        if k > j
            new_mask = max(new_mask, 0);
        end
    else
        new_mask = mask;
    end
end

function occ = get_occ(cur_im, tex_Pyr, j, cur_nnf, opt)
    occ = [];
    if isfield(opt, 'comp_penalty') && opt.comp_penalty > 0
        opt.scramble_nnf = 0;
        opt.comp_dist = 1;
        % todo use shorter input
        occ = distmex(cur_im(:, :, 1), tex_Pyr{j}(:, :, 1), cur_nnf, opt);
    end
end

function [mask] = get_mask(dist, occ, opt)
    bad = dist > 7;
    good = dist < 2;
    if nargin > 1 && ~isempty(occ)
        good = good & (occ < 2);
        bad = bad | (occ > 3);
    end
    % different mask modes
    if isfield(opt, 'bin_mask') && opt.bin_mask
        mask = single(-bad);
    else
        % fixed (perfect) >= 1
        % normal (good) in [0;1)
        % scramble (bad) < 0
        mask = single(good - bad);
    end
end
%% draw the angle
function [img] = drawAngleArrow(nnf)
    angleStep = 6;
    im_nnf = im2double(imread('test_bin.png'));
    img = rgb2lab(im_nnf);
    [w h] = size(nnf);
    img(:,:,1) = 0;
    for x = 1:w
        for y = 1:h
            index = nnf(y,x);
            img(y,x,2) = (index - angleStep)*127/angleStep;
            img(y,x,3) = (index - angleStep)*127/angleStep;
        end
    end
    img = lab2rgb(img);
end