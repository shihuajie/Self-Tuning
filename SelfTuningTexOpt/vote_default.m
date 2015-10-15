function [ out, weights ] = vote_default( nnf, src, P )
    if nargin < 3
        P = 10;
    end
    S = P - 1;
    szX = size(nnf, 1);
    szY = size(nnf, 2);
    nnfX = nnf(:, :, 1);
    nnfY = nnf(:, :, 2);
    num_ch = size(src, 3);
    % output of same class as source
    src_class = class(src);
    out = zeros(szX + S, szY + S, num_ch, src_class);
    weights = zeros(szX + S, szY + S, src_class);
    % dumb version
    for x = 1:szX
        for y = 1:szY
            px = nnfX(x, y) + 1;
            py = nnfY(x, y) + 1;
            bx = x:(x+S);
            by = y:(y+S);
            out(bx, by, :) = out(bx, by, :) + bilinearPatch( src, px:(px+S), py:(py+S) );
            weights(bx, by) = weights(bx, by) + 1;
        end
    end
    W = repmat(weights, [1, 1, num_ch]);
    mask = weights > 0;
    mask = repmat(mask, [1, 1, num_ch]);
    out(mask) = out(mask) ./ W(mask);
end

function [res] = bilinearPatch( src, px, py )
    px_f = floor(px);
    px_c = ceil(px);
    py_f = floor(py);
    py_c = ceil(py);
    tl = src(int32(px_f), int32(py_f), :);
    tr = src(int32(px_c), int32(py_f), :);
    br = src(int32(px_c), int32(py_c), :);
    bl = src(int32(px_f), int32(py_c), :);
    a = px(1) - px_f(1); a = 1 - a;
    b = py(1) - py_f(1); b = 1 - b;
    res = a * b * tl + (1-a) * b * tr + (1-a) * (1-b) * br + a * (1-b) * bl;
end