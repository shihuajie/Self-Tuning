function res = do_poisson_combination(im_gray, im_mask, grad_weight, im_dx, im_dy)

[h,w,nc] = size(im_gray);

im_dx(abs(im_dx) > 100) = 0;
im_dy(abs(im_dy) > 100) = 0;
ww = (1 - im_mask);
WG =  grad_weight * (ones(size(im_mask)) + 100 * ww);
WC =  ones(size(im_mask)) + 100 * ww;
[Gx,Gy] = getGMat(w,h,WG);
tmp = reshape(Gx * double(im_gray(:)), [h,w]);



[h,w,c] = size(im_gray);    

%mask = imerode(mask,ones(5));

n = size(Gx,1);
%disp(sprintf('using %f as prior in the grad stitched!', prior_weight));



Id = sparse( 1:n, 1:n, WC(:));
A = [Gx; Gy; Id];

b = [double(WG(:) .* im_dx(:)); double(WG(:) .* im_dy(:)); double(im_gray(:).* WC(:))];

% B = reshape(b, [h, w, 3]);
% imwrite(B(:, :, 1) / 10, 'gx.png');
% imwrite(B(:, :, 2) / 10, 'gy.png');
% imwrite(B(:, :, 3) / 100, 'L.png');

res = (A'*A) \ (A'*b);

res = reshape(res, [h,w]);

