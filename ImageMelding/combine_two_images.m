function [out_im] = combine_two_images(im1, im2, mask1, mask2, isMax)
nChannels = size(im1,3);
grad1 = mask1(:,:,1) .* (im1(:,:,4).^2 + im1(:,:,5).^2);
grad2 = mask2(:,:,1) .* (im2(:,:,4).^2 + im2(:,:,5).^2);

inds = find(repmat(grad1,[1 1 nChannels])  > repmat(grad2,[1 1 nChannels]) );
out_im = im2;
out_im(inds) = im1(inds);

%out_im(:,:,1) =  (mask1(:,:,1) .* im1(:,:,1) + mask2(:,:,1) .* im2(:,:,1))./(mask1(:,:,1)  + mask2(:,:,1));
%out_im(:,:,2) =  (mask1(:,:,1) .* im1(:,:,2) + mask2(:,:,1) .* im2(:,:,2))./(mask1(:,:,1)  + mask2(:,:,1));
%out_im(:,:,3) =  (mask1(:,:,1) .* im1(:,:,3) + mask2(:,:,1) .* im2(:,:,3))./(mask1(:,:,1)  + mask2(:,:,1));

%out_im(:,:,2:3) =  im1(:,:,2:3);
out_im(:,:,2) =  (mask1(:,:,1) .* im1(:,:,2) + mask2(:,:,1) .* im2(:,:,2))./(mask1(:,:,1)  + mask2(:,:,1));
out_im(:,:,3) =  (mask1(:,:,1) .* im1(:,:,3) + mask2(:,:,1) .* im2(:,:,3))./(mask1(:,:,1)  + mask2(:,:,1));
if isMax == 0
    out_im(:,:,4) =  (mask1(:,:,1) .* im1(:,:,4) + mask2(:,:,1) .* im2(:,:,4))./(mask1(:,:,1)  + mask2(:,:,1));
    out_im(:,:,5) =  (mask1(:,:,1) .* im1(:,:,5) + mask2(:,:,1) .* im2(:,:,5))./(mask1(:,:,1)  + mask2(:,:,1));
end