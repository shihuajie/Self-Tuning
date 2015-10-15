function [ F ] = imfill( I, color, modulo, thresh, shift )

    [szX, szY, ch] = size(I);
    [Y, X] = meshgrid(1:szY, 1:szX);
    P = find(mod(X + shift(1), modulo(1)) <= thresh | mod(Y + shift(2), modulo(2)) <= thresh);
    F = I;
    for c = 1:ch
        Fn = F(:, :, c);
        Fn(P) = color(c);
        F(:, :, c) = Fn;
    end
end

