function [ B, P, X, D ] = find_grid_base( I, varargin )
%FIND_GRID_BASE Find the best vectorial base to explain features in I

    % default options
    options.threshold = 0.5;
    options.min_dist = 5;
    options.min_hist = 5;
    options.plotting = 0;
    options.use_grid = 1;
    options.method = 'lattice';
    options.sigma = 3.0;
    options.penalty = 3.0;
    options.min_angle = pi/16;
    
    % user options
    if numel(varargin) == 1 && isstruct(varargin{1})
        % a big options variable
        new_opts = varargin{1};
        fields = fieldnames(new_opts);
        for f = 1:numel(fields)
            fname = fields{f};
            options.(fname) = new_opts.(fname);
        end
    elseif mod(numel(varargin), 2)
        % invalid!
        error('Invalid number of options parameters, use pairs!');
    else
        % pairs of elements
        for i = 1:2:numel(varargin)
            name = lower(varargin{i});
            value = varargin{i+1};
            options.(name) = value;
        end
    end
    
    % default values
    B = [0, 0; 0, 0];
    P = [0; 0];
    
    [szY, szX] = size(I);
    min_dist = options.min_dist;
    max_dist = [szX, szY] / 2;

    addpath(fullfile('..', 'vlfeat', 'toolbox', 'mex', mexext));
    addpath('Utils');

    % default
    B = [0, 0; 0, 0];
    P = [0; 0];

    [feat, desc] = vl_sift(I);
    N = size(desc, 2);
    scores = pdist(double(desc'));
    SP = squareform(scores < options.threshold * mean(scores)); % tril( . )
    [rows, cols] = find(SP);
    if numel(rows) < 3
      return
    end
    
    % get offsets
    offX = (feat(1, rows) - feat(1, cols)).';
    offY = (feat(2, rows) - feat(2, cols)).';
    sgnX = sign(offX); % offX = offX .* sgnX;
    sgnY = sign(offY); % offY = offY .* sgnY;
    R = ceil([max(offX), max(offY)]);
    
    % filter by distance to origin
    %   too close => noise data
    %   too far => meaningless
    valid = filter_by_distance(offX, offY, min_dist, max_dist);
    if numel(valid) < 1
      return
    end
    offX = offX(valid); offY = offY(valid);
    sgnX = sgnX(valid); sgnY = sgnY(valid);
    fx_ind = rows(valid); fy_ind = cols(valid);
    
    % filter using a 2d histogram
    if options.use_grid
        valid = filter_with_2dgrid( offX, offY, min_dist, options.min_hist );
        if numel(valid) < 1
          return
        end
        offX = offX(valid); offY = offY(valid);
        sgnX = sgnX(valid); sgnY = sgnY(valid);
        fx_ind = fx_ind(valid); fy_ind = fy_ind(valid);
    end
    
    if nargout > 2
        X = [offX, offY];
    end
    if nargout > 3
        D = [feat(1, :).' feat(2, :).'];
    end
    
    if isempty(offX) || isempty(offY)
        return;
    end

    if strcmp(options.method, 'lattice')
        addpath(fullfile('..', 'PatchMatch', 'dist'));
        if ~isfield(options, 'range')
            options.range = options.min_dist; 
        end
            
        [B, P] = latticemex( [offX, offY], options );
        B = [B(:, 2), B(:, 1)]; % x,y -> y,x
    else
        % find the offsets
        [bX, pX, DX] = find_base( offX, min_dist, options.method );
        [bY, pY, DY] = find_base( offY, min_dist, options.method );
        B = [0, 0];
        P = [0, 0];
        if isempty(bX) && isempty(bY)
            return
        end
        if ~isempty(bX)
            B(2) = bX;
            P(2) = pX;
        end
        if ~isempty(bY)
            B(1) = bY;
            P(1) = pY;
        end
    
        % plot stuff
        if options.plotting
            subplot(2, 1, 1); plot(DX(1, :), DX(2, :)); 
            title('X Offsets'); xlabel('Offset n'); ylabel('Normalized sigma');
            if pX, hold all; plot(bX, DX(2, bX - min_dist + 1), '*r'); hold off; end
            subplot(2, 1, 2); plot(DY(1, :), DY(2, :));
            title('Y Offsets'); xlabel('Offset n'); ylabel('Normalized sigma');
            if pY, hold all; plot(bY, DY(2, bY - min_dist + 1), '*r'); hold off; end
        end
    end
end

function [newX, newY] = rot2d(X, Y, theta)
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; % transformation matrix
    newXY = R * [X'; Y'];
    newX = newXY(1, :)';
    newY = newXY(2, :)';
end

function valid = filter_by_distance( offX, offY, min_dist, max_dist )
    valid = find(abs(offX) < max_dist(1) & abs(offY) < max_dist(2) ...
        & sqrt(offX .* offX + offY .* offY) > min_dist);
end

function valid = filter_with_2dgrid( offX, offY, w, min_count )
    % shift everything to first quadrant
    offX = offX - min(offX(:));
    offY = offY - min(offY(:));
    
    % divide grid into small parts
    gridX = 1 + w * round(offX / w);
    gridY = 1 + w * round(offY / w);
    
    % accumulate on the grid
    G = max([gridX, gridY])
    grid = accumarray([gridX(:), gridY(:)], 1.0, G);
    L = sub2ind(G, gridX, gridY);
    valid = find(grid(L) > min_count); % filtering here
end

function [b, p, D] = find_base( X, P0, method )
    if nargin < 2
        P0 = 5;
    end
    Pmax = round(max(X) / 1.7);
    S = zeros(1, Pmax - P0 + 1);
    K = zeros(1, Pmax - P0 + 1);
    for P = P0:Pmax
        c = P / 2;
        XmP = mod(X, P); % modulo
        XmP = XmP - (XmP >= c) * P; % centering
        
        % fit normal distribution
        npd = fitdist(XmP, 'normal');
        S(P - P0 + 1) = npd.sigma;
        K(P - P0 + 1) = kurtosis(XmP);
    end
    figure; multimshow(2, 2, @plot, S, S ./ (P0:Pmax), K);
    switch lower(method)
        case 'sigma'
            D = [P0:Pmax; S ./ (P0:Pmax)];
            [~, s_idx] = min(D(2, :));
            b = s_idx + P0 - 1;
            p = D(2, s_idx) <= 0.2;
        case 'kurtosis'
            D = [P0:Pmax; K];
            [~, k_idx] = max(D(2, :));
            b = k_idx + P0 - 1;
            p = D(2, k_idx) > 3;
        otherwise
            error('Unknown method %s', method);
    end
end
