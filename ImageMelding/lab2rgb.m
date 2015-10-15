
function [out_im] = lab2rgb(input_im)

C = makecform('lab2srgb');
input_im(isinf(input_im)) = 0;
out_im = applycform(double(input_im), C);
