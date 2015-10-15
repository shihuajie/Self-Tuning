%% hist_show Show the histogram change from a picture to another
%
%% Input
%   im_in     the source image (from)
%   im_out    the result image (to)
%   options   the texture generation options
%
%% Output
%   h         the set of histograms for the matched channels
%   f         the histogram figure handle (optional)
%
function [ h, f ] = create_hist( im_in, im_out, options )
    if(nargin == 1 && isstruct(im_in))
        h = im_in;
        N = numel(h.src);
        do_disp = 1;
    else
        N = numel(options.hist_channels);
        h.src = cell(N, 1);
        h.trg = cell(N, 1);
        h.bins = cell(N, 1);
        h.opt = options;
        h.src_sz = size(im_in);
        h.trg_sz = size(im_out);
        % bins
        for i = 1:N
            range = options.hist_ranges_max(i) - options.hist_ranges_min(i);
            binSize = range / options.hist_bins(i);
            minBin = options.hist_ranges_min(i) + binSize/2;
            maxBin = options.hist_ranges_max(i) - binSize/2;
            h.bins{i} = minBin:binSize:maxBin;
        end
        % histograms
        for i = 1:N
            imch = im_in(:, :, i); 
            h.src{i} = hist(imch(:), h.bins{i}); 
            imch = im_out(:, :, i); 
            h.trg{i} = hist(imch(:), h.bins{i}); 
        end
        do_disp = nargout == 2;
    end
    % display
    if do_disp
        f = figure('visible', 'off');
    else
        return
    end
    set(f, 'Color', 'w');
    for i = 1:N
        subplot(1, N, i);
        h1 = h.src{i} / sum(h.src{i});
        h2 = h.trg{i} / sum(h.trg{i});
        plot(h.bins{i}, h1, 'r', h.bins{i}, h2, 'b');
        axis square;
    end

end

