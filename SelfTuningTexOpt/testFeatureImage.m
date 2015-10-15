    %%modfied by huajie 2015-8-6
   % addpath(fullfile('..', 'PatchMatch', 'dist'));
   % addpath('Utils');

%     trg_gc_file = 'src_gc1.jpg'
%     trg_gc  = imread(trg_gc_file);
%     %trg_gc = grb2gray(trg_gc);
%     %trg_gc(:,:,2:3) = 0;
%     thresh = graythresh(trg_gc);     %自动确定二值化阈值；
%     trg_gc = im2bw(trg_gc,thresh);       %对图像自动二值化即可。
% 
%     %trg_gc = 1 - trg_gc;
%     %%trg_gc = improc(trg_gc, 'edge,n,otsu');
%     trg_gc = improc(trg_gc, 'rev,sbwdist');
%    trg_gc = improc(trg_gc, 'n');
%    % trg_gc = mat2gray(trg_gc);
%    trg_gc = im2double(trg_gc);
%    imshow(trg_gc);
%     
%     imwrite(trg_gc, 'trg_gc_feature.bmp','bmp');


%      fp =  fopen('weight.txt','r');
%      while ~feof(fp)     
%         %data = fscanf(fp,'%f')
%         
%         data = fgetl(fp)
%         d = str2num(data);
%         %params.weightIndex = data;
% %         try
% %             [res_im, res_NN] = synthesize(src_im, init_im, trg_mask, options);
% %         catch ex
% %             getReport(ex)
% %             % if it's not for debug, then we exit
% %             if ~isfield(options, 'debug')
% %               exit(1);
% %             end
% %         end
%      end
%      fclose(fp);


% image = imread('im.jpg');
% R = image(:,:,1);
% G = image(:,:,2);
% B = image(:,:,3);
% thresh = graythresh(R);
% R = im2bw(R,thresh);
% thresh = graythresh(G);
% G = im2bw(G,thresh);
% thresh = graythresh(B);
% B = im2bw(B,thresh);
% 
% color = 'cyan';
% isDilate = false;
% isBlur = true;
% 
% 
% switch color
%     case 'red'
%         out = R;
%     case 'yellow'
%         out = logical(R) & logical(G);
%     case 'green'
%         out = G;
%     case 'cyan'
%         out = logical(G) & logical(B);
%     case 'blue'
%         out = B;
%     case 'pink'
%         out = logical(R) & logical(B);
%     case 'white'
%         out = logical(R) & logical(G) & logical(B);
%     case 'black'
%         out = logical(R) & logical(G) & logical(B); 
% end
% 
% if(isDilate)
%     SE = strel('disk',5,0);
%     out = imdilate(out,SE);
% end
% if(isBlur)
%     out = improc(out,'gaussian');
% end
% 
% improc(out, 'n:0:100');
% imshow(out);
% 
% function out = myfunc(a)
% 
% out = a 
%     switch a
%         case 'red'
%             sprintf('%s\n','red');
%         case 'green'
%             sprintf('%s\n','green');
% % end
% im = imread('leaf.jpg');
% figure;
% subplot(2,8,1);
% imshow(im);
% im = im(:,end:-1:1,:);
% size(im)
% subplot(1,2,2);
%  imshow(im);
% figure;
% subplot(2,8,1);
% imshow(im);
%  title('原图');
%  for i= 0:7
%      B = imrotate(im, i*45, 'bilinear');
%      size(B)
%      subplot(2,8,i+1);
%      imshow(B);
%      f_im = B(:,end:-1:1,:);
%      subplot(2,8,i+9);
%      imshow(f_im);
%  end
% 
% angle=-45;                       %要旋转的角度，旋转方向为顺时针?
% img=imread('leaf.jpg');       %这里v为原图像的高度，u为原图像的宽度?
% imshow(img);                    %这里y为变换后图像的高度，x为变换后图像的宽度?
% [h w]=size(img);  
% theta=angle/180*pi; 
% rot=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];  
% pix1=[1 1 1]*rot;               %变换后图像左上点的坐标?
% pix2=[1 w 1]*rot;               %变换后图像右上点的坐标?
% pix3=[h 1 1]*rot;               %变换后图像左下点的坐标?
% pix4=[h w 1]*rot;               %变换后图像右下点的坐标? 
% height=round(max([abs(pix1(1)-pix4(1))+0.5 abs(pix2(1)-pix3(1))+0.5]));     %变换后图像的高度?
% width=round(max([abs(pix1(2)-pix4(2))+0.5 abs(pix2(2)-pix3(2))+0.5]));      %变换后图像的宽度?
% imgn=zeros(height,width);  
% delta_y=abs(min([pix1(1) pix2(1) pix3(1) pix4(1)]));            %取得y方向的负轴超出的偏移量?
% delta_x=abs(min([pix1(2) pix2(2) pix3(2) pix4(2)]));            %取得x方向的负轴超出的偏移量? 
% for i=1-delta_y:height-delta_y    
%     for j=1-delta_x:width-delta_x 
%         pix=[i j 1]/rot;                                %用变换后图像的点的坐标去寻找原图像点的坐标，                                                                       %否则有些变换后的图像的像素点无法完全填充?       ?
%         float_Y=pix(1)-floor(pix(1));        
%         float_X=pix(2)-floor(pix(2));             
%         if pix(1)>=1 && pix(2)>=1 && pix(1) <= h && pix(2) <= w                  
%             pix_up_left = [floor(pix(1)) floor(pix(2))];          %四个相邻的点?            
%             pix_up_right = [floor(pix(1)) ceil(pix(2))];             
%             pix_down_left = [ceil(pix(1)) floor(pix(2))];             
%             pix_down_right = [ceil(pix(1)) ceil(pix(2))];           
%             value_up_left = (1-float_X)*(1-float_Y);              %计算临近四个点的权重?            
%             value_up_right = float_X*(1-float_Y);             
%             value_down_left = (1-float_X)*float_Y;             
%             value_down_right = float_X*float_Y;                                              
%             imgn(i+delta_y,j+delta_x) = value_up_left*img(pix_up_left(1),pix_up_left(2))+ ... 
%                                         value_up_right*img(pix_up_right(1),pix_up_right(2))+ ...                                        ?
%                                         value_down_left*img(pix_down_left(1),pix_down_left(2))+ ...                                        ?
%                                         value_down_right*img(pix_down_right(1),pix_down_right(2));
%         end
%     end
% end
% figure,imshow(uint8(imgn))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
r = 125;
img=im2double(imread('leaf.jpg'));
im = img;
[h w]=size(im(:,:,1)); 
width = ceil(((h*h +w*w)^(1/2)));
src_im = zeros(width, width, 3);
h_edge = floor((width-h)/2);
w_edge = floor((width-w)/2);
src_im(h_edge+1:h_edge+h, w_edge+1:w_edge+w, :) = im;
imshow(src_im);
r = width/2;

new_im = zeros(width,width,3);
figure;
for index = 0:15
      theta=mod(index, 8) * pi/4; 
      rot=[ cos(theta)  sin(theta) ;
           -sin(theta)  cos(theta) ];
      for dx = 1:width
          for dy = 1:width
              
            px = dx - r;
            py = dy - r;
            if(index >= 8)
                px(1) = px(1)*(-1);
            end

            tx = rot*[px py]';
            nx = r + tx(1);
            ny = r + tx(2);
            if (nx < 1 || nx > width || ny < 1 || ny > width) continue;end
%             new_im(dx,dy,:)  = src_im(round(nx), round(ny),:);
            pix = [nx ny];
            float_Y=pix(1)-floor(pix(1));        
            float_X=pix(2)-floor(pix(2));
            pix_up_left = [floor(pix(1)) floor(pix(2))];       
            pix_up_right = [floor(pix(1)) ceil(pix(2))];             
            pix_down_left = [ceil(pix(1)) floor(pix(2))];             
            pix_down_right = [ceil(pix(1)) ceil(pix(2))];           
            value_up_left = (1-float_X)*(1-float_Y); 
            value_up_right = float_X*(1-float_Y);             
            value_down_left = (1-float_X)*float_Y;             
            value_down_right = float_X*float_Y;                                              
            new_im(dx,dy,:)  = value_up_left*src_im(pix_up_left(1),pix_up_left(2),:)+ ... 
                                        value_up_right*src_im(pix_up_right(1),pix_up_right(2),:)+ ...                                        ?
                                        value_down_left*src_im(pix_down_left(1),pix_down_left(2),:)+ ...                                        ?
                                        value_down_right*src_im(pix_down_right(1),pix_down_right(2),:);
          end
      end
      subplot(4, 4, index+1)
      imshow(im2double(new_im));
end
