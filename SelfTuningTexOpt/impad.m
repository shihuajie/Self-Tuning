function padded = impad(img, P_h, P_w)
    [h, w, c] = size(img);
    if nargin < 3
        P_w = P_h;
    end
    if islogical(img)
        rev_log = 1;
        img = single(img);
    else
        rev_log = 0;
    end
    padded = zeros(h + P_h, w + P_w, c, 'like', img);
    padded(1:(end - P_h), 1:(end - P_w), :) = img;
    % padding
    if P_h > 0
        padded((end - P_h + 1):end, 1:(end - P_w), :) = repmat(img(end, :, :), P_h, 1);
    end
    if P_w > 0
        padded(1:(end - P_h), (end - P_w + 1):end, :) = repmat(img(:, end, :), 1, P_w);
    end
    if P_h > 0 && P_w > 0
        padded((end - P_h + 1):end, (end - P_w + 1):end, :) = repmat(img(end, end, :), P_h, P_w);
    end
    % type transform
    if rev_log
        padded = logical(padded);
    end
end