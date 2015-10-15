function output_im = make_blend_SShape(input_im, gamma)
%input_im = alpha_mask;
%gamma = 2;
L = 10000;
T = drawS(L + 1, gamma);
transforms = int32(round((L -1) * input_im)) + 1;
%output_im = input_im;

output_im = T(transforms);