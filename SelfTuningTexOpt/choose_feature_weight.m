function flag = choose_feature_weight( I )
    if ischar(I)
        I = im2double(imread(I));
    end
    
    addpath(fullfile('..', 'PatchMatch', 'dist'));
    addpath('Utils');
    
    timer = tic;
    
    % by default, use high weight
    flag = 1;

    % get Lab channels
    Lab = improc(I, 'lab');
    
    % 3D color accumulation using gaussian kernels
    channels = [100, 0, 100, 100, -80, 80, 100, -80, 80];
    opt.sigma = 3;
    E = improc(I, 'edge,n,otsu');
    F0 = improc(E, 'rev,sbwdist');
    F = improc(F0, 'n:0:100');
    Edist = F(E);
    zeroDist = ceil(max(Edist(:)));
    maxDist = max(zeroDist + 7, min(2 * zeroDist, zeroDist + 20));
    for c = 1:3:7
        opt.channels = channels;
        opt.channels(c) = 1;
        opt.mask = F > zeroDist + 3 & F <= maxDist;
        H = reshape(colorhistmex(Lab, opt), [100, 100]);
        
        % check the number of different colors
        % if the area is too big, then it's a dangerous case
        % since even though one color might be predominant, there are many
        % other colors to consider
        area = numel(find(H > 0));
        if area > 5000
            flag = -0.5;
            break;
        end
        
        % threshold the low counts
        maxH = max(H(:));
        B = H > maxH / opt.sigma;
        H1 = H .* B;
        % [~, N] = bwlabel(B);
        
        % count number of local maxima
        [x, y, ~, t] = imextrema(H1);
        N = numel(find(t == 1));
        
        % if extrema are too close, it doesn't count
        D = max(pdist([x(t == 1), y(t == 1)]));
        if D < 15
            continue;
        end
        
        % multiple => we should not use a high weight
        %   since this would imply mixing!
        if N > 1
            flag = -N;
            break;
        end
    end
    
    duration = toc(timer);
    
    save_data(duration, 'feature_elapsed.dat');
end
