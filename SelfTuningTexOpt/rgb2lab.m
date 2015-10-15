function [ out_im ] = rgb2lab( in_im )
    C = makecform('srgb2lab');
    out_im = applycform(in_im, C);
end

