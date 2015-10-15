function [H W C T] = imsize(img)

H = size(img,1);
W = size(img,2);
C = size(img,3);
T = size(img,4);

if nargout==1
	H = [H W C T];
end

	
	