%FIND_ALIGNMENT - find the best translation of a picture B onto A
%
%% Input
%   A           the fixed image
%   B           the image to best translate
%   maskA       the mask of valid pixels in A
%   maskB       the mask of valid pixels in B
%   validity    a binary function for valid positions in C
%   fast        whether to speed up processing
%
%% Output
%   tX      the best translation of B to maximize correlation with A
%   tY
%   maxK    the maximum correlation that was found
%
function [ tX, tY, maxK ] = find_alignment( A, B, maskA, maskB, validity, fast )

    % not 1 until it is tested
    if nargin < 6
        fast = 0;
    end

    [szAX, szAY, nCh] = size(A);
    [szBX, szBY, ~] = size(B);
    
    % relative shift
    relShX = 0;
    relShY = 0;

    % default arguments (and simplification)
    if nargin < 3 || isempty(maskA)
        maskA = ones(szAX, szAY) == 1;
    else
        % cut the bottom and right parts that are useless
        [szAX, szAY, A, maskA, relShX, relShY] = cut_useless(A, maskA, fast);
    end
    if nargin < 4 || isempty(maskB)
        maskB = ones(szBX, szBY) == 1;
    else
        [szBX, szBY, B, maskB] = cut_useless(B, maskB);
    end
    szCX = szAX + szBX - 1;
    szCY = szAY + szBY - 1;
    if nargin < 5 || isempty(validity)
        valid = ones(szCX, szCY) == 1;
    else
        [X, Y] = ndgrid(1:szCX, 1:szCY);
        valid = logical(bsxfun(@(tx, ty)validity(tx - szBX + relShX, ty - szBY + relShY), X, Y));
    end

    C = zeros(szCX, szCY);
    % per channel
    for c = 1:nCh
        % channels
        chA = A(:, :, c);
        chB = B(:, :, c);
        % mean value
        chMean = mean([chA(maskA); chB(maskB)]);
        % mean-centering
        chA = chA - chMean;
        chA(~maskA) = 0;
        chB = chB - chMean; % note: use the same mean!
        chB(~maskB) = 0;

        % cross-correlation
        C = C + xcorr2(chA, chB);
    end
    K = C .* valid;
    [maxK, ind] = max(K(:));
    [shX, shY] = ind2sub(size(C), ind(1));
    tX = shX - szBX + relShX;
    tY = shY - szBY + relShY;
    1;
end

function [szAX, szAY, newA, newMaskA, shAX, shAY] = cut_useless(A, maskA, cutBefore)
    if nargin < 3
        cutBefore = 0;
    end
    [szAX, szAY] = size(maskA);
    [Y, X] = meshgrid(1:szAY, 1:szAX);
    offAX = X .* maskA;
    offAY = Y .* maskA;
    maxAX = max(offAX(:));
    maxAY = max(offAY(:));
    
    % cutting the left / top to optimize speed
    if cutBefore
        minAX = min(offAX(offAX > 0));
        minAY = min(offAY(offAY > 0));
    else
        minAX = 1;
        minAY = 1;
    end
    % relative offset
    shAX = minAX - 1;
    shAY = minAY - 1;
    
    % cutting helps speed up xcorr2
    newA = A(minAX:maxAX, minAY:maxAY, :);
    newMaskA = maskA(minAX:maxAX, minAY:maxAY);
    szAX = maxAX - minAX + 1;
    szAY = maxAY - minAY + 1;
end
