function [ res_im, res_NN ] = synth_func( file, out_dir, varargin )
%SYNTH_FUNC Start a full texture synthesis
%
% SYNOPSIS
%   [ im, nnf ] = synth_func( file, out_dir, ... )
%
% INPUT
%   - file      exemplar filename
%   - out_dir   directory to output results in
%   - ... parameters for the synthesis
%
% OUTPUT
%   - res_im    the result image
%   - res_NN    the result nnf
%

    addpath(fullfile('..', 'PatchMatch', 'dist'));
    addpath(fullfile('..', 'ImageMelding'));
    addpath('Utils');

    %% Global parameters
    global out_directory; %directory name for saving intermediate files
    if isempty(out_directory)
        [~, name, ~] = fileparts(file);
        if isempty(out_dir)
            out_directory = name;
        else
            out_directory = [out_dir filesep name];
        end
        mkdir(out_directory);
    elseif ~exist(out_directory, 'dir')
        mkdir(out_directory);
    end
    
    %% Parameters
    options = get_synth_options(file, varargin{:})

    %% Call string
    % serializing the call to store it
    call = [ 'synth_func( ''' file ''', ''' out_dir ''' ' ];
    save_call = 1;
    if ~isempty(varargin)
        argCount = size(varargin, 2);
        for i = 1:argCount
            arg = varargin{i};
            if isstruct(arg)
                save_call = 0;
            end
            call = [ call ', ' toMString(arg) ];
        end
    end
    call = [ call ' )' ];
    if save_call
        save_data(call, 'call');
    end
    
    %% Texture exemplar (source)
    src_im = im2double(imread(file)); % load source texture
    imwrite(src_im, fullfile(out_directory, 'exemplar.png'));
    if size(src_im, 3) == 1
      src_im = repmat(src_im, [1, 1, 3]); % gray to rgb
      options.nnf_channels = [1]; % no need to use a or b
    end
    src_im = build_image( src_im, options, 'src' );
    
    %% Texture target
    if isfield(options, 'target')
        if isnumeric(options.target)
            trg_im = create_target(options.target);
        elseif ischar(options.target)
            % load target file
            trg_im = im2double(imread( options.target ));
        else
            error('Invalid target data');
        end
    else
        if isfield(options, 'size')
            trg_im = create_target(options.size);
        else
            trg_im = create_target(options.targetSize);
        end
    end
    trg_im = build_image( trg_im, options, 'trg' );

    %% Target mask
    trg_mask = [];
    if isfield(options, 'mask')
        mask_str = options.mask;
        mask_mode = '+';
        if ~exist(mask_str, 'file')
            mask_mode = mask_str(1);
            if ~exist(mask_str(2:end), 'file')
                error('Mask "%s" is not a valid file.', mask_mode);
            end
            mask_str = mask_str(2:end);
        end
        trg_mask = im2double(rgb2gray(imread(mask_str)));
        switch mask_mode
          case '-'
            trg_mask = -single(trg_mask); % 1 => should be randomized
            trg_fill = 1; % the rest should be fixed
          case '+'
            trg_mask = single(trg_mask); % 1 => should be fixed
            trg_fill = -1; % the rest should be randomized
          otherwise
            error('Invalid mask mode "%s"', mask_mode);
        end
        % fill the rest using the distance transform
        if ~isfield(options, 'bin_mask') || ~options.bin_mask
            mask_margin = options.patch_size; % by default the patch size
            if isfield(options, 'mask_margin')
                mask_margin = options.mask_margin;
            end
            trg_mask = trg_mask + trg_fill * single( bwdist(trg_mask * trg_fill < 0) > mask_margin );
        end
    end

    %% Class is ... single!
    trg_im = single(trg_im);
    src_im = single(src_im);
    
    %% Randomization
    rng_seed = 'shuffle'; % default to use the timestamp as seed
    if isfield(options, 'rand_seed')
       rng_seed = options.rand_seed;
    end
    rng_gen = 'twister'; % default to use Mersenne Twister
    if isfield(options, 'rand_generator')
        rng_gen = options.rand_generator;
    end
    rng(rng_seed, rng_gen);
    
    %% Lattice generators
    if ischar(options.smart_init) || options.smart_init || options.aligned_search
        if ischar(options.smart_init)
            % load cached result
            G = load_mat(options.smart_init, 'G');
            options.generators = G;
        else
            % compute lattice generators
            grid_opts.measure = 'fmeasure';
            grid_opts.penalty = 3;
            grid_opts.approximate = 1; % be reasonable, the exact version is overkill
            grid_opts.beta = 0.1;
            grid_opts.use_grid = 1; % get as many points as possible
            grid_opts.range = [10, 400]; % reasonable generators

            [G, P] = find_grid_base(src_im(:, :, 1), grid_opts)
            options.generators = [G, P];

            % if nothing, then use patch size
            if numel(find(G)) == 0
                G = options.patch_size * eye(2, 2);
                options.generators = [G, P];

                % disable aligned search, as it's meaningless without
                % generators
                options.aligned_search = 0;
            end
        end
    end
    
    %% modfied by huajie 2015-8-6
    src_gc_file = 'src_gc.png';
    trg_gc_file = 'trg_gc.png';
    if ~isempty(out_directory)
        src_gc_file = strcat(out_directory, '\', src_gc_file);
        trg_gc_file = strcat(out_directory, '\', trg_gc_file);
    end
    
    src_gc = single(imread(src_gc_file));
    trg_gc = single(imread(trg_gc_file));
    src_lab = rgb2label(src_gc);
    trg_lab = rgb2label(trg_gc);
      
    src_im = cat(3, src_im, src_lab);
    trg_im = cat(3, trg_im, trg_lab);
    %% end
    
    %% Initial random image
    if options.smart_init
        if ischar(options.smart_init)
            % load cached result
            data = load_mat(options.smart_init, 1);
            init_im = data.init_im;
            G = data.G;
            options.generators = G;
            t = 0;
        else
            % initialize the target with smart blocks
            t = tic;
            init_im = smart_init(trg_im, src_im, options.generators(1:2, 1:2), options);    
            init_im = single(init_im);
            t = toc(t);
            fprintf('Smart init in %f seconds.\n', t);

            if isfield(options, 'smart_init_cache') && options.smart_init_cache
                G = options.generators;
                save(fullfile(out_directory, 'smart_init.mat'), ...
                    'init_im', 'G');
            end
        end
        
        % save result for posterity
        I0 = lab2rgb(init_im(:, :, 1:3));
        if isfield(options, 'smart_init_only') && options.smart_init_only
            imwrite(I0, fullfile(out_directory, 'img.png')); % the final result
            save_data(t, 'elapsed');
            save_data('done');
            res_im = I0;
            res_NN = [];
            return
        else
            imwrite(I0, fullfile(out_directory, 'smart_init.png'));
            save_data(t, 'smart_init_elapsed');
        end
        
    elseif ~isfield(options, 'random_init') || options.random_init
        % random blocks
        if isfield(options, 'rand_block')
            rb = options.rand_block;
            P_x = rb(1);
            if numel(rb) > 1
                P_y = rb(2);
            else
                P_y = P_x;
            end
            int_loc = 1;
        else
            P_x = options.patch_size;
            P_y = options.patch_size;
            int_loc = 0;
        end
        % initial image by selecting random patches from the exemplar
        init_im = random_init(trg_im, src_im, P_x, P_y, int_loc);
        init_im = single(init_im);
        
        % be dumb?
        if isfield(options, 'dumb_init') && options.dumb_init
            range = min([size(trg_im); size(src_im)], [], 1);
            init_im(1:range(1), 1:range(2), :) = src_im(1:range(1), 1:range(2), :);
        end
        
    else
        init_im = trg_im;
    end
    
    % keep the fixed parts from the original target if there is a mask
    if ~isempty(trg_mask)
        num_ch = size(trg_im);
        mask = repmat(trg_mask >= 1, [1, 1, num_ch]);
        init_im(mask) = trg_im(mask);
    end
    % keep the fixed channels ... fixed
    if options.fixed_channels
        init_im(:, :, options.fixed_channels) = trg_im(:, :, options.fixed_channels);
    end
    
%% added by huajie 2015-9-11
%% start
%      src_gc = imread(src_gc_file);
%     thresh = graythresh(src_gc);     %自动确定二值化阈值；
%     src_gc = im2bw(src_gc,thresh);       %对图像自动二值化即可。
%     src_gc = 1 - src_gc;
%     src_gc = improc(src_gc, 'rev,sbwdist');
%     src_gc = improc(src_gc, 'n');
%     imwrite(src_gc, fullfile(out_directory, 'src_gc_feature.bmp'));
%     src_im = cat(3, src_im, src_gc);
%     
%      trg_gc = imread(trg_gc_file);
%     thresh = graythresh(trg_gc);     %自动确定二值化阈值；
%     trg_gc = im2bw(trg_gc,thresh);       %对图像自动二值化即可。
%     trg_gc = 1 - trg_gc;
%     trg_gc = improc(trg_gc, 'rev,sbwdist');
%     trg_gc = improc(trg_gc, 'n');
%     imwrite(trg_gc, fullfile(out_directory, 'trg_gc_feature.bmp'));
%     init_im = cat(3, init_im, trg_gc);

%     src_gc = im2double(imread(src_gc_file));
%     src_gc_R = src_gc(:,:,1);
%     SE = strel('disk',5,0);
%     src_gc_R = imdilate(src_gc_R,SE);
%     src_gc_R = improc(src_gc_R,'gaussian');
%     src_gc(:,:,1) = src_gc_R;
%     src_gc = improc(src_gc,'n:0:100');
%     src_im = cat(3, src_im, src_gc);
%  
%     trg_gc = im2double(imread(trg_gc_file));
%     trg_gc_R = trg_gc(:,:,1);
%     trg_gc_G = trg_gc(:,:,2);
%     SE = strel('disk',5,0);
%     trg_gc_R = imdilate(trg_gc_R,SE);
%     trg_gc_G = imdilate(trg_gc_G,SE);
%     trg_gc_R = improc(trg_gc_R,'gaussian');
%     trg_gc_G = improc(trg_gc_G,'gaussian');
%     trg_gc(:,:,1) = trg_gc_R;
%     trg_gc(:,:,2) = trg_gc_G;
%     trg_gc = improc(trg_gc,'n:0:100');
%     init_im = cat(3, init_im, trg_gc);

% %     src_gc = imread(src_gc_file);
% %     src_gc_R = getColorChannel(src_gc, 'red', true, true);
% %     src_gc(:, :, 1) = src_gc_R;
% %     src_im = cat(3, src_im, src_gc);
% %     
% %     trg_gc = imread(trg_gc_file);
% %     trg_gc_R = getColorChannel(trg_gc, 'red', true, true);
% %     trg_gc(:, :, 1) = trg_gc_R;
% %     trg_gc_C = getColorChannel(trg_gc, 'cyan', false, false);
% %     trg_gc(:, :, 2) = trg_gc_C;
% %     init_im = cat(3, init_im, trg_gc);
     
%       src_gc_image = mat2gray(src_lab);
%       imwrite(src_gc_image, [out_directory filesep 'src_lab.png']);
%       trg_gc_image = mat2gray(trg_lab);
%       imwrite(trg_gc_image, [out_directory filesep 'trg_lab.png']);
%% end
%% Synthesis
    %% added by huajie 2015-9-11
    % process more weight
    start_time = tic;
    fp =  fopen([out_directory filesep 'weight.txt'],'r');
    while ~feof(fp)     
        data = str2num(fgetl(fp));
        outPath = strcat('constraint_weight_',num2str(data));
        out_directory = [out_dir filesep name filesep outPath];
        if ~exist(out_directory, 'dir')
            mkdir(out_directory);
        end
        options.weight = data;
        try
            [res_im, res_NN] = synthesize(src_im, init_im, trg_mask, options);
        catch ex
            getReport(ex)
            % if it's not for debug, then we exit
            if ~isfield(options, 'debug')
              exit(1);
            end
        end
    end
    fclose(fp);
    % elapsed time
    total_elapsed_time = toc(start_time)
    out_directory = [out_dir filesep name];
    save_data(total_elapsed_time, 'total_elapsed_time');
%% end
%     try
%         [res_im, res_NN] = synthesize(src_im, init_im, trg_mask, options);
%     catch ex
%         getReport(ex)
%         % if it's not for debug, then we exit
%         if ~isfield(options, 'debug')
%           exit(1);
%         end
%         
%     end
end

function [res] = create_target( target )
    switch numel(target)
      case 1
        res = zeros(target, target, 3);
      case 2
        res = zeros(target(1), target(2), 3);
      otherwise
        res = target;
    end
end

function I = load_im( fname )
    if isnumeric(fname)
        % it's an actual image! no need to do anything
        I = fname;
    elseif strcmp(fname(end-2:end), 'mat')
        % load mat file
        I = load_mat(fname);
    else
        % load image
        I = im2double(imread(fname));
    end
end

function [res] = build_image( im, options, prefix )
    res = rgb2lab(im); % srgb to Lab
    szX = size(im, 1);
    szY = size(im, 2);
    
    % we add the Gx, Gy channels
    res = cat(3, res, zeros(szX, szY, 2));
    
    % gradient computation
    [res_gx, res_gy] = gradient1(res(:, :, 1), options.patch_grad_weight);
    res(:, :, options.grad_channels) = cat(3, res_gx, res_gy);
    
    % dimensions
    [H, W, ~] = size(res);
    
    % look for extra channels
    i = 1;
    while isfield(options, ['src_extra' num2str(i)])
        if isfield(options, [prefix '_extra' num2str(i)])
            extra = load_im( options.([prefix '_extra' num2str(i)]) );
            extra = improc(extra, [prefix '_extra' num2str(i)], options);
        else
            extra = zeros(szX, szY);
        end
        
        % channel weight
        if isfield(options, ['weight_extra' num2str(i)])
            w = options.(['weight_extra' num2str(i)]);
            extra = extra * w;
        end
        
        % checking the dimensions
        [h, w] = size(extra);
        if h ~= H || w ~= W
            error('Invalid sizes: %s is %d by %d, whereas extra%d is %d by %d', prefix, H, W, i, h, w); 
        end
        
        % add channel at the end
        res = cat(3, res, extra);
        
        i = i + 1;
    end
end

function [str] = toMString( data )
    if isnumeric(data)
        if numel(data) < 2
            str = num2str(data);
        else
            if ~isa(data, 'double')
                str = [ class(data) '(['];
            else
                str = '[';
            end
            num_rows = size(data, 1);
            num_cols = size(data, 2);
            for r = 1:num_rows
                str = [ str num2str(data(r, 1)) ];
                for c = 2:num_cols
                   str = [ str ', ' num2str(data(r, c)) ];
                end
                if r < num_rows
                   str = [ str '; '];
                end
            end
            str = [ str ']'];
            if ~isa(data, 'double')
                str = [ str ')'];
            end
        end
    elseif ischar(data)
        str = ['''' data ''''];
    elseif isstruct(data)
        str = 'struct(...)';
    else
        fprintf('Unsupported data type for serialization: %s\n', class(data));
        str = '__unsupported__';
    end
end

%% get one channel about a color, also have options to dilate or blured
function out = getColorChannel(image, color, isDilate, isBlur)
    R = image(:,:,1);
    G = image(:,:,2);
    B = image(:,:,3);
    thresh = graythresh(R);
    R = im2bw(R,thresh);
    thresh = graythresh(G);
    G = im2bw(G,thresh);
    thresh = graythresh(B);
    B = im2bw(B,thresh);
    
    switch color
        case 'red'
            out = R;
        case 'yellow'
            out = logical(R) & logical(G);
        case 'green'
            out = G;
        case 'cyan'
            out = logical(G) & logical(B);
        case 'blue'
            out = B;
        case 'pink'
            out = logical(R) & logical(B);
        case 'white'
            out = logical(R) & logical(G) & logical(B);
        case 'black'
            out = logical(R) & logical(G) & logical(B); 
    end
    out = single(out)
    if(isDilate)
        SE = strel('disk',5,0);
        out = imdilate(out,SE);
    end
    if(isBlur)
        out = improc(out,'gaussian');
    end
    improc(out, 'n:0:100');
    imshow(out);
end

%% added by huajie 2015-9-15
% rgb image transform to a single channel label to represent different
% annotations
function [ out] = rgb2label(image)
%     image = double(imread(imagepath));
    [gcSzX, gcSzY, ~] = size(image);
    out = zeros(gcSzX, gcSzY,'single');
    label = round(image/128);
    for x = 1:gcSzX
        for y = 1:gcSzY
            t = label(x, y, :);
            out(x, y) = [t(1) t(2) t(3)]*[9 3 1]';
        end
    end
end
