function [ gradX, gradY ] = gradient1( img, weight )
%GRADIENT1 Compute a 1 channel simple difference gradient

    szX = size(img, 1);
    szY = size(img, 2);
    numCh = size(img, 3);
    if numCh ~= 1
        error('Expect a 1 channel image!');
    end
    
    [Gx, Gy] = getGMat(szY, szX, ones(szX, szY));
    % [src_gx, src_gy] = gradient(src_im(:, :, 1)); <-- not the good one!
    I = double(reshape(img, szX * szY, 1));
    gradX = reshape(Gx * I, szX, szY);
    gradY = reshape(Gy * I, szX, szY);
    
    % use weight
    if nargin >= 2 && weight > 0 && weight ~= 1
        gradX = gradX * weight;
        gradY = gradY * weight;
    end

end

