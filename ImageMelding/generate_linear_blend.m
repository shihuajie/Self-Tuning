% function [alpha_mask] = generate_linear_blend(input_mask)
% 
% x = input_mask(1,:,:);
% st_ind = find(x, 1, 'first');
% end_ind = find(x, 1, 'last');
% y = ones(size(x));
% l = end_ind - st_ind;
% y(st_ind: end_ind - 1) = 1 - (1:l)/ l;
% y(end_ind: end) = 0;
% 
% alpha_mask = repmat(y, [size(input_mask,1),1]);


function [alpha_mask] = generate_linear_blend(input_mask)

x = input_mask(1,:,:);
st_ind = find(x, 1, 'first');
end_ind = find(x, 1, 'last');
y = ones(size(x));
l = end_ind - st_ind;
M = drawS(l, 1);
y(st_ind: st_ind +length(M) - 1) = 1 - M;%(1:l)/ l;
y(st_ind +length(M): end) = 0;
N = 100;

alpha_mask = repmat(y, [size(input_mask,1),1]);


