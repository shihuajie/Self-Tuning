function [ I ] = improc( I0, key, varargin )
%IMPROC Process an image with filters
%
% Input
%   - I0        the image to process with filters
%   - key       the options key prefix
%   - options   the options struct
%
% Output
%   - I         the processed image
%   
    I = I0; % by default

    p = {};
    if numel(varargin) > 0 && isstruct(varargin{1})
        options = varargin{1};
        % load the proc_str data
        if isfield(options, [key '_proc'])
            proc_str = options.([key '_proc']);
        else
            proc_str = 'gray';
        end

        % get parameters
        p_i = 1;
        while isfield(options, sprintf('%s_proc_p%d', key, p_i))
            p{p_i} = options.(sprintf('%s_proc_p%d', key, p_i));
            p_i = p_i + 1;
        end
    else
        proc_str = key;
        p = varargin;
    end
    % substitute parameters
    if numel(p) > 0
        proc_str = sprintf(proc_str, p{:});
    end

    % split processing string
    % strsplit(proc_str, ',') is only available from Matlab 2013a
    proc = strparts(proc_str, ',');
    % apply each filter
    for i = 1:numel(proc)
        parts = strparts(proc{i}, ':');
        switch(lower(parts{1}))
            case {'f', 'gray', 'flatten'}
                % flatten to gray (if needed)
                if size(I, 3) > 1
                    I = rgb2gray(I);
                end
                
            case {'unflatten', 'ungray'}
                I = repmat(I, [1, 1, 3]);
                
            case 'lab'
                % rgb to Lab
                I = rgb2lab(I);
                
            case 'rgb'
                % Lab to rgb
                I = lab2rgb(I);
                
            case {'n', 'norm', 'normalize'}
                % normalize:a:b
                if numel(parts) < 2
                    a = 0;
                else
                    a = str2double(parts{2});
                end
                if numel(parts) < 3
                    b = 1;
                else
                    b = str2double(parts{3});
                end
                I_min = min(I(:));
                I_max = max(I(:));
                I = a + (I - I_min) / (I_max - I_min) * (b - a);
                
            case {'h', 'hist', 'histnorm'}
                % histogram normalization
                if numel(parts) < 2
                    thresh = 0.01;
                else
                    thresh = str2double(parts{2});
                end
                I = imadjust(I, stretchlim(I, thresh), []);
                
            case {'r', 'rev', 'reverse'}
                % image reverse
                I = 1 - I;
                
            case '*'
                % multiplication
                if numel(parts) < 2
                    factor = 1.0;
                else
                    factor = str2double(parts{2});
                end
                I = I * factor;
                
            case '+'
                % increment
                if numel(parts) < 2
                    incr = 0.0;
                else
                    incr = str2double(parts{2});
                end
                I = I + incr;
                
            case {'du', 'downup'}
                if numel(parts) < 2
                    factor = 0.5;
                else
                    factor = 1/str2double(parts{2});
                end
                if numel(parts) < 3
                    method = 'lanczos3';
                else
                    method = parts{3};
                end
                
                sz = [size(I, 1), size(I, 2)];
                I = imresize(imresize(I, ceil(sz * factor), method), sz, method);
                
            case {'g', 'normal', 'gaussian'}
                % g:sigma:kSize
                if numel(parts) < 2
                    sigma = 1;
                else
                    sigma = str2double(parts{2});
                end
                if numel(parts) < 3
                    kSize = 3;
                else
                    kSize = str2double(parts{3});
                end
                
                I = imfilter(I, fspecial('gaussian', kSize, sigma), 'replicate');
                
            case {'e', 'edge', 'edgemap'}
                % edgemap:nms
                addpath(fullfile('..', 'EdgeMap'));
                addpath(fullfile('..', 'PiotrToolbox', 'channels'));
                
                % TODO implement!!!
                model = load_mat(fullfile('..', 'EdgeMap', 'models', 'forest', 'modelFinal.mat'));
                model.opts.nPos=5e5; 
                model.opts.nNeg=5e5;
                model.opts.useParfor=0;
                model.opts.multiscale=0;
                model.opts.nTreesEval=4;
                model.opts.nThreads=4;
                if numel(parts) > 1
                    model.opts.nms = str2double(parts{2});
                else
                    model.opts.nms = 0;
                end
                I = edgesDetect(I, model);
                
            case {'d', 'grad', 'gradient'}
                [Dx, Dy] = gradient1(I);
                I = sqrt(Dx .* Dx + Dy .* Dy);
                
            case 'ggdt'
                % gddt:mu:delta:useI0:num_it
                addpath(fullfile('..', 'PatchMatch', 'dist'));
                
                if numel(parts) > 1
                    mu = str2double(parts{2});
                else
                    mu = 70;
                end
                if numel(parts) > 2
                    gamma = str2double(parts{3});
                else
                    gamma = 0;
                end
                if numel(parts) > 3
                    useI0 = strcmp(parts{4}, '1');
                else
                    useI0 = 0;
                end
                if numel(parts) > 4
                    num_it = str2double(parts{5});
                else
                    num_it = 5;
                end
                
                if useI0
                    I = ggdtmex(I, mu, I0, gamma, num_it);
                else
                    I = ggdtmex(I, mu, I, gamma, num_it);
                end
                
            case 'sggdt'
                % sgddt:mu:delta:useI0:num_it
                addpath(fullfile('..', 'PatchMatch', 'dist'));
                
                if numel(parts) > 1
                    mu = str2double(parts{2});
                else
                    mu = 70;
                end
                if numel(parts) > 2
                    gamma = str2double(parts{3});
                else
                    gamma = 0;
                end
                if numel(parts) > 3
                    useI0 = strcmp(parts{4}, '1');
                else
                    useI0 = 0;
                end
                if numel(parts) > 4
                    num_it = str2double(parts{5});
                else
                    num_it = 5;
                end
                
                if useI0
                    I = ggdtmex(I, mu, I0, gamma, num_it) - ggdtmex(1 - I, mu, 1 - I0, gamma, num_it);
                else
                    I = ggdtmex(I, mu, I, gamma, num_it) - ggdtmex(max(I(:)) - I, mu, 1 - I, gamma, num_it);
                end
                
            case 'bwdist'
                I = bwdist(~I); % bwdist takes from 1, not 0!
                
            case 'sbwdist'
                I = bwdist(~I) - bwdist(I);
                
            case 'otsu'
                if numel(parts) > 1
                    delta = str2double(parts{2});
                else
                    delta = 0;
                end
                thresh = graythresh(I) + delta;
                if numel(parts) > 2
                    min_thresh = str2double(parts{3});
                    thresh = max(min_thresh, thresh);
                end
                if numel(parts) > 3
                    max_thresh = str2double(parts{4});
                    thresh = min(max_thresh, thresh);
                end
                I = im2bw(I, thresh);
                
            case {'t-', 'thresh-', 'threshold-'}
                %threshold-:max
                thresh = str2double(parts{2});
                I = I .* (I <= thresh);
                 
            case {'t+', 'thresh+', 'threshold+'}
                thresh = str2double(parts{2});
                I = I .* (I >= thresh);
                 
            case '<'
                thresh = str2double(parts{2});
                I = double(I < thresh);
                 
            case '<='
                thresh = str2double(parts{2});
                I = double(I <= thresh);
                
            case '>'
                thresh = str2double(parts{2});
                I = double(I > thresh);
                
            case '>='
                thresh = str2double(parts{2});
                I = double(I >= thresh);
                
            case 'slic'
                %slic:regularizer:regionSize
                addpath(fullfile('..', 'vlfeat', 'toolbox', 'mex', mexext));
                if numel(parts) > 1
                    regularizer = str2double(parts{2});
                else
                    regularizer = 1;
                end
                if numel(parts) > 2
                    regionSize = str2double(parts{3});
                else
                    regionSize = 10;
                end
                
                [H, W, ch] = size(I);
                
                segments = vl_slic(single(I), regionSize, regularizer) + 1;
                S = max(segments(:));
                segment_groups = accumarray(segments(:), 1:numel(segments), [S, 1], @(x){x});
                Ich = reshape(I, H * W, ch);
                for g = 1:S
                    group_mean = mean(Ich([segment_groups{g}], :), 1);
                    for c = 1:ch
                        Ich([segment_groups{g}], c) = group_mean(c);
                    end
                end
                I = reshape(Ich, [H, W, ch]);
                
            case 'upscale'
                targetSize = str2double(parts{2});
                [H, W, ~] = size(I);
                scale = targetSize / min([H, W]);
                if numel(parts) > 2
                    method = parts{3};
                else
                    method = 'lanczos3';
                end
                I = imresize(I, scale, method);
                
            case 'rescale'
                if numel(parts) > 2
                    method = parts{3};
                else
                    method = 'lanczos3';
                end
                I = imresize(I, [size(I0, 1), size(I0, 2)], method);
        end
    end

end

function parts = strparts(str, delim)
    if exist('strsplit', 'file') > 0
        parts = strtrim(strsplit(str, delim));
    else
        parts = strtrim(strread(str, '%s', 'delimiter', delim));
    end
end