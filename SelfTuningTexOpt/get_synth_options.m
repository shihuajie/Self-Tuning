%% get_synth_options - process synthesis options
%
%% Input
%   set of arguments pairs param, value
%
%% Output
%
%   options   the full default options, with overwritten parameters
%             optimized when parameters can simplify other parameters
function [ options ] = get_synth_options( file, varargin )

    %% Default parameters
    % patch and distance types
    options.patch_size = 10; % patch size, area=patch_size * patch_size
    %modified by huajie
    options.patch_type = 'basicf'; % basic, basicf or affine
    % options.multires = [zeros(1, 2), ones(1, 5), ones(1, 3) * 2];
    % options.threads_num = 1;
    options.dist_type = 'ssd';
    
    % affine bounds
    options.min_scale = 0.9;% minimum scale in GPM 
    options.max_scale = 1.2;% maximum scale in GPM
    options.min_angle = -pi/2;%minimum rotation
    options.max_angle = pi/2;%maximum rotation

    % iteration parameters
    options.iterations = [10 * ones(5, 1), 7 * ones(5, 1)]; %[(42-(1:4)), (30-2*(5:10))]; % number of GPM (Generalized PatchMatch) iterations
    options.it_start_order = 0;
    options.rand_search = 6;
    options.rand_search_until = 4;
    % options.it_until_conv = 1;
    % options.coh_rand_search = 1;
    options.max_rand_search = 5;
    % options.mixed_rand_search = 1;
    % options.window_size = 0;
    % options.window_size_decr = 1;
    
    % completeness enforcement
    % options.comp_penalty = 100;
    % options.comp_weight = 2;
    options.comp_exponent = 1;
    % options.comp_prop = 0;
    options.comp_prop_until = 4;
    options.incomp_search_until = 4;
    % options.comp_abs_threshold = 0;
    % options.comp_rel_threshold = 1;
    
    % pyramid and scales
    options.num_scales = 10;
    options.num_em = 30; % number of EM iteration at coarsest scale
    options.decrease_factor = 3; % decrease of num_em per scale
    options.min_em = 5; % minimum number of EM
    options.coarsest_scale = 35; % minimum side for the pyramids
    options.is_log_upscale = 1;
    
    % automatic methods
    options.adaptive = 0;
    
    % bias and gain
    % options.min_Bias = single([-10   ; -50   ;   -50;      0; 0; 0]);%minimum bias in Lab space.  The last two channels are related to gradients.  As gradients can be drived from colors, if you leave it as it's the algorithm automatically uses gain of color channels.
    % options.max_Bias = single([10 ; 50 ;  50; 0;0; 0]);%maximum bias 
    % options.min_Gain = single([0.9  ; 0.9     ;   0.9;      1;1; 1]);%minimum gain
    % options.max_Gain = single([1.5  ; 1.5   ;   1.5;    1;1; 1]);%maximum gain

    % weights
    options.grad_weight = 0;
    options.feat_weight = 75;
    options.patch_grad_weight = 0; % the weight of the gradient in the patch distance

    % channels
    options.nnf_channels = 1:3; % [1:3, 6] % which channels to use for the nnf search: L,a,b,F
    options.feature_channels = []; % 6;
    options.grad_channels = 4:5;
    options.fixed_channels = []; % 7
    
    %% Smart initialization
    options.smart_init = 0;
    options.aligned_search = 0;
    options.max_aligned_search = 5;
    options.aligned_search_until = 4;

    %% Output parameters
    options.saveEach = 0; % whether to save everything at each scale
    options.saveChannels = 0; % save voted channels individually
    options.saveConvData = 0;
    options.saveDist = 0;
    options.saveEnv = 0;
    options.saveGradient = 0;
    options.saveHist = 0;
    options.saveLastNNF = 0;
    options.saveMask = 0;
    options.saveNNC = 1;
    options.saveNNF = 0;
    options.saveOccRatio = 0;

    %% Voting options
    % options.vote_weight = single( 1000 * fspecial('gaussian', P, 3));
    % options.vote_weight = single(ones(options.patch_size, options.patch_size));
    options.vote_method = 'default';
    % options.vote_channels = 1:5; % L,a,b,Gx,Gy
    % meanshift
    options.meanshift_window = single([20, 40, 40, 40, 40]);
    % histograms
    options.hist_bins = single([100, 20, 20]);
    options.hist_channels = single([0, 1, 2]);
    options.hist_ranges_min = single([0, -80, -80]);
    options.hist_ranges_max = single([100, 80, 80]);
    options.hist_weights = single([1000, 10, 10]);
    
    %% Texture parameters
    % options.interpolation = 'nn'; % 'bilinear' by default
    options.boundaries = 'same';
    options.weight          = 0; % added by huajie
    %% Overridden parameters
    if numel(varargin) > 1
        argCount = size(varargin, 2);
        if mod(argCount, 2) == 1
            error('Invalid number of arguments: %d should be even!', argCount);
        end
        for i = 1:2:argCount
            argName = varargin{i};
            argValue = varargin{i + 1};
            options.(argName) = argValue;
        end
    elseif numel(varargin) == 1
        % user-provided basic options
        user_options = varargin{1};
        if ~isstruct(options)
            error('Second option argument to get_synth_options is not a struct!');
        end
        % overwrite default with these
        user_fields = fieldnames(user_options);
        for i = 1:numel(user_fields)
            f = user_fields{i};
            options.(f) = user_options.(f);
        end
    end
    
    %% Shortcuts
    if options.adaptive
        options.extra_channels = 2;
        options.hist_params = 'default';
    end
    nnf_weights = ones(1, numel(options.nnf_channels));
    for i = 1:numel(options.nnf_channels)
        c = options.nnf_channels(i);
        opt_name = ['nnf_weight' num2str(c)];
        if isfield(options, opt_name)
            nnf_weights(i) = options.(opt_name);
            options.nnf_weights = nnf_weights;
        end
    end
    % automatic weighting
    if isfield(options, 'nnf_weights') && strcmp(options.nnf_weights, 'auto')
        if within_array(6, options.nnf_channels)
            F_flag = choose_feature_weight(file);
            dir = 1;
            if isfield(options, 'nnf_weights_reverse') && options.nnf_weights_reverse
                dir = -1;
            end
            if F_flag * dir >= 0
                Fweight = 3;
                if isfield(options, 'high_weight')
                    Fweight = options.high_weight;
                end
                nnf_weights(options.nnf_channels == 6) = Fweight;
            end
        end
         options.nnf_weights = nnf_weights; 
        %%options.nnf_weights = [1, 1, 1, 0, 1, 1, 1];%%added by huajie 2015-9-5
    end
    if isfield(options, 'extra_channels')
        [~, name, ext] = fileparts(file);
        edgeFile = fullfile('Results', 'edges', [name ext]);
        if exist(edgeFile, 'file')
            from = 'gray';
        else
            edgeFile = file;
            from = 'edge';
        end
        % F - c#6 - signed distance
        if options.extra_channels > 0
            options.src_extra1 = edgeFile;
            options.src_extra1_proc = [from ',n,otsu,rev,sbwdist,n:0:100'];
        end
        % E - c#7 - binary edge map
        if options.extra_channels > 1
            options.src_extra2 = edgeFile;
            options.src_extra2_proc = [from ',n,otsu,n:0:100'];
        end
    end
    if isfield(options, 'hist_params') && strcmp(options.hist_params, 'default')
        hist_params = [];
        for c = options.vote_channels
            switch c
                case 1
                    hist_params = [hist_params; 1 100 0 100 1];
                case 2
                    hist_params = [hist_params; 2 100 -80 80 1];
                case 3
                    hist_params = [hist_params; 3 100 -80 80 1];
                otherwise
                    if isfield(options, 'binary_channels') && within_array(c, options.binary_channels)
                        hist_params = [hist_params; c 2 0 100 1];
                    else
                        hist_params = [hist_params; c 100 0 100 1];
                    end
            end
        end
        hist_params
        options.hist_params = hist_params;
    end
    
    %% Interpretation
    F = fieldnames(options);
    for i = 1:numel(F)
        f = F{i};
        if ischar(options.(f)) && ~isempty(strfind(options.(f), '%'))
            options.(f) = interpret(options.(f), file, options);
        end
    end
    
    %% Transformations
    if isfield(options, 'vote_weight')
        options.vote_weight = parseWeights(options.vote_weight, options.patch_size);
    end
    if isfield(options, 'vote_weight_then')
        options.vote_weight_then = parseWeights(options.vote_weight_then, options.patch_size);
    end
    
    %% Simplifications
    if ~isfield(options, 'vote_channels')
        if options.grad_weight <= 0
            options.vote_channels = 1:3; % no need to vote for the gradient
        else
            options.vote_channels = 1:5; % vote the gradient
        end
    end

end

function res = interpret( str, file, options )
    [~, filename, ext] = fileparts(file);
    res = strrep( str, '%FILENAME', [filename ext] );
    res = strrep( res, '%FILEPART', filename );
    res = strrep( res, '%FILE', file );
    if isfield(options, 'target') && ischar(options.target)
        res = strrep( res, '%TARGET', options.target);
    end
    % replace options
    m = regexp( res, '%\{[a-zA-Z_\}', 'match' );
    for i = 1:numel(m)
        m_str = m{i};
        opt_name = m_str(3:end-1);
        if isfield( options, opt_name )
            opt_val = options.(opt_name);
            if ischar(opt_val)
                rep_str = opt_val;
            elseif isnumeric(opt_val)
                rep_str = num2str(opt_val);
            else
                continue;
            end
            res = strrep( res, m_str, rep_str );
        end
    end
end

function [weights] = parseWeights(value, P)
    switch class(value)
        case 'char'
            value = lower(value);
            if(value(1) == 'n') % centered gaussian
                sigma = str2double(value(2:end));
                weights = single(fspecial('gaussian', P, sigma) * 100); % small values are bad
            elseif(value(1) == 't') % top-left gaussian
                sigma = str2double(value(2:end));
                weights = single(fspecial('gaussian', 2 * P, sigma) * 100);
                weights = weights(P+1:end, P+1:end);
            elseif(value(1) == 'd') % double gaussian (TL + BR)
                sigma = str2double(value(2:end));
                weights = single(fspecial('gaussian', 2 * P, sigma) * 100);
                weights = weights(P+1:end, P+1:end) + weights(1:P, 1:P);
            elseif(value(1) == 'q') % quadruple gaussian (each corner)
                sigma = str2double(value(2:end));
                weights = single(fspecial('gaussian', 2 * P, sigma) * 100);
                weights = weights(P+1:end, P+1:end) + weights(1:P, 1:P) ...
                    + weights(P+1:end, 1:P) + weights(1:P, P+1:end);
            else
                error('Unsupported vote_weight: %s', value);
            end

        case 'single'
            weights = value;
            
        otherwise
            weights = single(value);
    end
end

function b = within_array(elem, arr)
    b = ~isempty(find(arr(:) == elem, 1));
end
