function [out_im] = do_init_guess(in_im, in_mask)

  in_mask_d = imdilate(in_mask,strel('square',3));

out_im = in_im;
num_channels = size(in_im,3);
for ch = 1 : num_channels
    tmp = in_im(:,:,ch);
%     for kk = 0 : 10
%        
%         rL = tmp(1 + kk,:);
%         rR = tmp(end - kk,:);
%         rT = tmp(:,1 + kk);
%         rB = tmp(:,end - kk);
%         
%         rL_M = in_mask(1 + kk,:);
%         rR_M = in_mask(end - kk,:);
%         rT_M = in_mask(:,1 + kk);
%         rB_M = in_mask(:,end - kk);
%         
%         
%         
%         
%         rL = interp1(find(rL_M == 0),rL(find(rL_M == 0)),[1:length(rL)],'pchip');
%         rR = interp1(find(rR_M == 0),rR(find(rR_M == 0)),[1:length(rR)],'pchip');
%         rT = interp1(find(rT_M == 0),rT(find(rT_M == 0)),[1:length(rT)],'pchip');
%         rB = interp1(find(rB_M == 0),rB(find(rB_M == 0)),[1:length(rB)],'pchip');
%         
%         tmp(1 + kk,:) = rL;
%         tmp(end - kk ,:) = rR;
%         tmp(:,1 + kk) = rT;
%         tmp(:,end - kk) = rB;
%         
%         in_mask(1 + kk,:) = zeros(size(rL_M));
%         in_mask(end - kk,:) = zeros(size(rR_M));
%         in_mask(:,1 + kk) = zeros(size(rT_M));
%         in_mask(:,end - kk) = zeros(size(rB_M));
%     end
    out_im(:,:,ch) = roifill(tmp, in_mask_d);
    
    
    
    

    
    
end
