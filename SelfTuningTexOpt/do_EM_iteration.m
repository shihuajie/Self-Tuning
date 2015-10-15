function [output_im, output_NNF, output_dist, output_iter, output_conv] = do_EM_iteration(...
    trg_im, tex_im_PYR, org_im, ...
    start_iter, num_em, scale_ind, ...
    nnf, mask, options)
%DO_EM_ITERATION - run one EM iteration at a given scale of PM + Poisson
%
% Input
%   - trg_im      the target synthesized image
%   - tex_im_PYR  the texture source pyramid
%   - org_im      the original image (synthesis)
%   - start_iter  id of the starting iteration
%   - num_em      the number of em to do
%   - scale_ind   the scale index
%   - nff         the previous NNF to use
%   - mask        the image mask to use (fixed: >= 1, scramble: < 0)
%   - options     the options
%
% Output
%   - output_im   the new synthesized image
%   - output_NNF  the new set of nnfs
%   - output_dist the result patch distance
%   - output_iter the next iteration id
%   - output_conv the convergence data
%
  
global out_directory;

addpath(fullfile('..', 'PatchMatch', 'dist'));
addpath('Utils');

%% Sizes
num_ch = size(org_im, 3);
tex_im = tex_im_PYR{scale_ind, 1}; % the current texture

%% Dist feature specific
if within_array(6, options.vote_channels)
    F0 = tex_im(:, :, 6);
    F0 = improc(F0, 'n:0:100');
    tex_im(:, :, 6) = F0;
end

%% Input
if options.multires > 0
  tex_input = cell(1, 1 + options.multires);
  for r = 1:(1 + options.multires)
      tex_input{r} = tex_im_PYR{scale_ind, r}(:, :, options.nnf_channels);
  end
else
  tex_input = tex_im(:, :, options.nnf_channels);
end

%% NNF Normalization
if isfield(options, 'nnf_norm_channels')
    % set tex_input
    for c = 1:numel(options.nnf_norm_channels)
        channel = tex_input(:, :, c);
        if c == 2 || c == 3
            c_min = -80;
            c_max = 80;
        else
            c_min = 0;
            c_max = 100;
        end
        this_min = min(channel(:));
        this_max = max(channel(:));
        tex_input(:, :, c) = c_min + (channel - this_min) * ( ...
            (c_max - c_min) / (this_max - this_min) ...
        );
    end
end

%% NNF Weighting
nnf_weights = [];
if isfield(options, 'nnf_weights')
    nnf_weights = options.nnf_weights;
    if numel(nnf_weights) == 1
        nnf_weights = repmat(nnf_weights, 1, numel(options.nnf_channels));
    elseif numel(nnf_weights) ~= numel(options.nnf_channels)
        error('The size of nnf_weights does not match the size of nnf_channels!');    
    end
    % set tex_input
    for c = 1:numel(options.nnf_channels)
        tex_input(:, :, c) = nnf_weights(c) * tex_input(:, :, c);
    end
end
    
%% NNF Mask for voting
if ~isempty(mask)
    % nnf mask
    P = options.patch_size;
    nnf_mask = mask(1:end - P + 1, 1:end - P + 1);
    nnf_mask(nnf_mask >= 1) = 0.5; % do not fix the nnf (let propagation happen)
    options.vote_method = 'weighted';
    options.vote_mask = single(bwdist(nnf_mask >= 0));
    
    % linear mask for the image reset (not NNF!)
    L_mask = mask >= 1;
    trg_mask = repmat(L_mask, [1, 1, num_ch]);
else
    nnf_mask = [];
    L_mask = [];
    trg_mask = [];
end

%% EM iterations

% output the initial 0 iteration
if options.saveEach
    if options.saveChannels
        for c = options.vote_channels
            imwrite(uint8(trg_im(:, :, c)) + 1, jet(101), fullfile(out_directory, ...
                sprintf('init-%03d_c%d_s%02d.png', start_iter, c, scale_ind)));
        end
    else
        imwrite(lab2rgb(trg_im(:,:,1:3)), [out_directory filesep ...
            sprintf('init-%03d_s%02d.png', start_iter, scale_ind) ]);
    end
end

% convergence data
field_count = 5;
field_list = { 'prop_count', 'rs_count', 'as_count', 'is_count', 'max_dist', 'mean_dist', 'coh_ratio' };
if options.saveOccRatio
    field_list{6} = 'occ_ratio';
end
if (nargout >= 3)
    for i = 1:field_count
        output_conv.(field_list{i}) = [];
    end
end

% bidir similarity
if iscell(nnf)
    rev_nnf = nnf{2};
    nnf = nnf{1};
else
    rev_nnf = -1;
end

% do at most num_em EM steps
iter = start_iter;
confidence = [];
for e = 1 : num_em
    % special flag
    last_em = e == num_em;

    % reset special channels
    if options.patch_grad_weight > 0
        % we recompute the gradient
        [trg_gx, trg_gy] = gradient1(trg_im(:, :, 1), options.patch_grad_weight);
        trg_im(:, :, options.grad_channels) = cat(3, trg_gx, trg_gy);
    end
    if options.fixed_channels
        trg_im(:, :, options.fixed_channels) = org_im(:, :, options.fixed_channels);
    end

    %% Feature specific
    if within_array(6, options.vote_channels)
        F = trg_im(:, :, 6);
        minF = min(F(:));
        if minF > 0
            trg_im(:, :, 6) = F - minF;
        end
    end

    %% NNF Update
    opt = options;

    % augmenting the nnf
    if isfield(opt, 'ext_nnf')
        if opt.ext_nnf && opt.max_rand_search > 0
            extOpt = rmfield(opt, 'ext_nnf');
            extOpt.max_rand_search = 5;
            extOpt.self_nnf = 1;
            extOpt.min_patch_disp = options.patch_size;
            % extOpt.iterations = 20;
            extNNF = nnmex(tex_input, tex_input, [], [], extOpt);

            % saveNNF(extNNF(:, :, 1:2), size(tex_im), iter, scale_ind, e, 'ext-');
            % saveNNC(extNNF(:, :, 1:2), size(tex_im), iter, scale_ind, e, 'ext-');

            % put it in the options for the real nnf
            opt.ext_nnf = extNNF;
        else
            opt = rmfield(opt, 'ext_nnf');
        end
    end

    % using a weighted ssd
    if strcmp(opt.dist_type, 'wssd')
        opt.weight_type = 'variance';
        if ~isempty(nnf)
            P = opt.patch_size;
            % create distance weight
            if isempty(confidence)
                b = segmentmex( nnf, opt.patch_type );
                b = impad(single(b), P - 1);
                b = imshift([P, P] - round(P / 2), b);
                w = 1 - imfilter(b, fspecial('gaussian', ceil(P / 2), P / 2), 'replicate');
            else
                % confidence is here the pixel variance
                % confidence = 0 => sure
                % confidence > 0 => unsure
                w = 1 ./ (1 + confidence);
            end
            w = single(w);
            % normalize to [0;1]
            range = max(w(:)) - min(w(:));
            w = (w - min(w(:))) / range;
            % shift to [0.001;1]
            min_w = 0.1;
            range = 1.0 - min_w;
            w = (w * range) + min_w;

            % w: confidence in the pixel
            % w=1 where it's good
            % w<1 where it's bad
            opt.dist_weight = w;
            % imagesc(w); colorbar;
            % drawnow;
        else
            % not yet possible since every pixel is equally likely
            opt.dist_type = 'ssd';
        end
    end

    % iteration-dependent parameters
    if isfield(opt, 'it_until_conv') && opt.it_until_conv >= 1 && e < num_em
        opt.it_until_conv = 0; % only for the last one
    end

    if isfield(opt, 'it_start_order')
        opt.it_start_order = 1 - (opt.it_start_order == 1);
    end
    if isfield(options, 'rand_seed')
        opt.rand_seed = randi(intmax('int32'), 'int32');
    end

    trg_channels = trg_im(:, :, opt.nnf_channels);
    if opt.multires > 0
        trg_input = create_multires(trg_channels, opt.multires, opt.patch_size, opt.multires_weight);
    else
        trg_input = trg_channels;
    end
    if isfield(options, 'nnf_norm_channels')
        % set tex_input
        for c = 1:numel(options.nnf_norm_channels)
            channel = trg_input(:, :, c);
            if c == 2 || c == 3
                c_min = -80;
                c_max = 80;
            else
                c_min = 0;
                c_max = 100;
            end
            this_min = min(channel(:));
            this_max = max(channel(:));
            trg_input(:, :, c) = c_min + (channel - this_min) * ( ...
                (c_max - c_min) / (this_max - this_min) ...
            );
        end
    end
    if ~isempty(nnf_weights)
        for c = 1:numel(options.nnf_channels)
            trg_input(:, :, c) = nnf_weights(c) * trg_input(:, :, c);
        end
    end

    % reverse nnf
    if with_reverse_nnf(options)
        % different options as the reverse exemplar (the target) is
        % changing until the end of the synthesis, so we need more search
        rev_opt = opt;
        rev_opt.rand_search = 10; % always search
        rev_opt.coh_rand_search = 1; % search even when coherent as things change quickly
        rev_opt.comp_penalty = 0; % no completeness penalty here!
        rev_nnf = nnmex(tex_input, trg_input, rev_nnf, [], rev_opt);
        opt.reverse_nnf = rev_nnf;
        if options.saveNNF && options.saveEach
            saveNNF(rev_nnf(:, :, 1:2), size(trg_im), iter, scale_ind, e, 'rev-');
        end
        if options.saveNNC && options.saveEach
            saveNNC(rev_nnf(:, :, 1:2), size(trg_im), iter, scale_ind, e, 'rev-');
        end
    end

    %% EM iteration
    [newNNF, dist, cur_vote] = em_iteration({trg_input, trg_im}, {tex_input, tex_im}, nnf, nnf_mask, opt);

    % stop randomization
    if isfield(options, 'scramble_nnf')
        options = rmfield(options, 'scramble_nnf');
    end
    % possibly output the distance
    if options.saveDist && options.saveEach
        saveDist(dist, iter, scale_ind, e);
    end
    % possibly output the nnf map
    if options.saveNNF && options.saveEach
        saveNNF(newNNF(:, :, 1:2), size(tex_im), iter, scale_ind, e);
    end
    % possibly output texture patch use count
    if options.saveNNC && options.saveEach
        saveNNC(newNNF(:, :, 1:2), size(tex_im), iter, scale_ind, e);
    end


    %% Updating the data
    trg_im(:, :, opt.vote_channels) = cur_vote(:, :, :);
    nnf = newNNF;

    if isfield(options, 'dF') && options.dF
        fprintf('Fixing the F channel.\n');
        F = trg_im(:, :, 6);
        optF = options;
        if isfield(options, 'dF_rand_search')
            optF.rand_search = options.dF_rand_search;
        end
        optF.hist_params = [1, 100, 0, 100, 1];
        optF.nnf_channels = 1;
        optF.vote_channels = 1;
        [~, ~, F] = em_iteration({F, F}, {F0, F0}, newNNF, [], optF);
        trg_im(:, :, 6) = F;
    end

    % force fixed region
    if ~isempty(mask)
        trg_im(trg_mask) = org_im(trg_mask);
    end

    % debugging the vote
    if isfield(options, 'debug') && options.debug > 1
        voted = vote_default(newNNF, tex_im(:, :, 1:3), opt.patch_size);
        imwrite(lab2rgb(voted), [out_directory filesep ...
          sprintf('img-%03d_s%02d_e%02d-def.png', iter, scale_ind, e) ]);

        if options.saveGradient
            gxi = options.grad_channels(1);
            gyi = options.grad_channels(2);
            trg_gx = trg_im(:, :, gxi);
            trg_gy = trg_im(:, :, gyi);
            save([out_directory filesep ...
              sprintf('Gx-%03d_s%02d_e%02d.mat', iter, scale_ind, e) ], 'trg_gx');
            save([out_directory filesep ...
              sprintf('Gy-%03d_s%02d_e%02d.mat', iter, scale_ind, e) ], 'trg_gy');
        end
    end

    %% Gradient integration
    if options.grad_weight
        % the channels of interest
        gxi = options.grad_channels(1);
        gyi = options.grad_channels(2);

        % Poisson equation solving
        res = do_poisson_combination(trg_im(:,:,1), ones(size(trg_im(:,:,1))), ...
          options.grad_weight, trg_im(:,:,gxi), trg_im(:,:,gyi));
        trg_im(:, :, 1) = res;

        % forcing the fixed parts
        if ~isempty(mask)
            trg_im(L_mask) = org_im(L_mask);
        end

        % Clamping brightness
        % if scale_ind >= 9
        %    res(res > max_val_B) = max_val_B;
        %    trg_im(:,:,1) = res;
        % else
        %    res(res > 0.95 * max_val_B) = max_val_B;
        %    trg_im(:,:,1) = res;
        % end
    end

    % possibly save the histograms
    if options.saveHist && (e == 1 || last_em)
        histData = create_hist(tex_im, trg_im, options);
        save(fullfile(out_directory, sprintf('hist-%03d_s%02d_e%02d.mat', iter, scale_ind, e) ), ...
            'histData');
        % export_fig(sprintf(['%s' filesep 'hist%02d_Iter%02d.png'], ...
        %    out_directory, scale_ind, e), '-m2');
    end

    % output the current result
    if options.saveEach
        if options.saveChannels
            for c = options.vote_channels
                imwrite(uint8(trg_im(:, :, c)) + 1, jet(101), fullfile(out_directory, ...
                    sprintf('img-%03d_c%d_s%02d_e%02d.png', start_iter, c, scale_ind, e)));
            end
        else
            imwrite(lab2rgb(trg_im(:, :, 1:3)), fullfile(out_directory, ...
                sprintf('img-%03d_s%02d_e%02d.png', iter, scale_ind, e) ));
        end
    end

    % we went through a new iteration
    iter = iter + 1;

end

if options.saveChannels
    opts = options;
    opts.vote_channels = 1:3;
    opts.vote_method = 'default';
    lab = votemex(trg_im(:, :, 1:3), tex_im(:, :, 1:3), newNNF, opts);
    imwrite(lab2rgb(lab), fullfile(out_directory, sprintf('img-lab_s%02d.png', scale_ind) ));
end

% saving the last if needed
if ~options.saveEach && options.last_scale
    % intermediate data
    if options.saveDist,  saveDist(dist), end
    if options.saveNNF,   saveNNF(nnf, size(tex_im)); end
    if options.saveNNC,   saveNNC(nnf, size(tex_im)); end
    % result image
    if options.saveChannels
        for c = options.vote_channels
            imwrite(uint8(trg_im(:, :, c)) + 1, jet(101), fullfile(out_directory, ...
                sprintf('img_c%d.png', c)));
        end
    else
        imwrite(lab2rgb(trg_im(:, :, 1:3)), fullfile(out_directory, 'img.png'));
    end
end
% saving convergence data for the last scale
if options.last_scale
    % convergence data
    save_data(max(dist(:)), 'max_dist');
    save_data(mean(dist(:)), 'mean_dist');
    
    % environment save
    if options.saveEnv
        save([out_directory filesep sprintf('em_env_s%02d_i%03d.mat', ...
                scale_ind, iter)]);
    end
end

output_im = trg_im;
output_NNF = nnf;
output_dist = dist;
output_iter = iter;

if with_reverse_nnf(options)
    output_NNF = {nnf, rev_nnf};
end

end

function [newNNF, dist, vote] = em_iteration(trg, tex, nnf, nnf_mask, options)
    trg_input = trg{1};
    tex_input = tex{1};
    trg_im = trg{2};
    tex_im = tex{2};
    
    %% NNF Update
    
    if isfield(options, 'rand_seed')
        options.rand_seed = randi(intmax('int32'), 'int32');
    end
    
    % NNF Search
    % [newNNF, dist, conv] = nnmex(trg_input, tex_input, nnf, nnf_mask, options)
    [newNNF, dist] = nnmex(trg_input, tex_input, nnf, nnf_mask, options);

    %% NNF Voting

    if isfield(options, 'rand_seed')
        options.rand_seed = randi(intmax('int32'), 'int32');
    end

    % NNF Vote
    fprintf('Voting ...\n');
    %if isfield(options, 'weight_type') && strcmp(options.weight_type, 'variance')
    %    [vote, conf] = votemex(trg_im(:, :, options.vote_channels), ...
    %                   tex_im(:, :, options.vote_channels), newNNF, options);
    %else
    vote = votemex(trg_im(:, :, options.vote_channels), ...
                       tex_im(:, :, options.vote_channels), newNNF, options);
                
    fprintf('... done!\n');
end

function b = with_reverse_nnf(options)
    if isfield(options, 'reverse_nnf')
        rev_nnf = options.reverse_nnf;
        if numel(rev_nnf) > 1
            b = 1;
        else
            b = rev_nnf ~= 0;
        end
    else
        b = 0;
    end
end

function [] = saveDist(dist, iter, scale_ind, e, prefix)
    if nargin < 5
        prefix = '';
    end
    global out_directory;
    img = dist / max(dist(:));
    if nargin == 1
        imwrite(img, fullfile(out_directory, 'dist.png'));
    else
        imwrite(img, [out_directory filesep ...
          sprintf('dist-%s%03d_s%02d_e%02d.png', prefix, iter, scale_ind, e) ]);
    end
end

function [] = saveNNF(nnf, tex_sz, iter, scale_ind, e, prefix)
    if nargin < 6
        prefix = '';
    end
    global out_directory;
    tex_szX = tex_sz(1);
    tex_szY = tex_sz(2);
    % normalize the nnf
    nnf(:, :, 1) = nnf(:, :, 1) / tex_szX;
    nnf(:, :, 2) = nnf(:, :, 2) / tex_szY;
    img = cat(3, nnf(:,:,1:2), zeros(size(nnf, 1), size(nnf, 2)));
    if nargin == 2
        imwrite(img, fullfile(out_directory, 'nnf.png'));
    else
        imwrite(img, [out_directory filesep ...
          sprintf('nnf-%s%03d_s%02d_e%02d.png', prefix, iter, scale_ind, e) ]);  
    end
end

function [] = saveNNC(nnf, tex_sz, iter, scale_ind, e, prefix)
    if nargin < 6
        prefix = '';
    end
    global out_directory;
    tex_szX = tex_sz(1);
    tex_szY = tex_sz(2);
    % accumulate data
    nnX = 1 + floor(nnf(:, :, 1));
    nnY = 1 + floor(nnf(:, :, 2));
    if numel(find(nnX <= 0)) > 0 || numel(find(nnY <= 0)) > 0
       fprintf('Invalid nnX / nnY'); 
    end
    nnC = accumarray([nnX(:), nnY(:)], 1.0, [tex_szX, tex_szY]);
    img = nnC / max(nnC(:));
    if nargin == 2
        imwrite(img, fullfile(out_directory, 'nnc.png'));
    else
        imwrite(img, [out_directory filesep ...
          sprintf('nnc-%s%03d_s%02d_e%02d.png', prefix, iter, scale_ind, e) ]);
    end
end

function b = within_array(elem, arr)
    b = ~isempty(find(arr(:) == elem, 1));
end
