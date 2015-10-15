function [h] = im(x, varargin)
% IM - Display an image
%
%	[h] = im(x)
%	[h] = im(x, clip)
%
%
%IN:
%	x - HxWxC the image to be displayed
%	clip - if is present, clip the image data to [0..1] instead of
%		scaling. If it is a 2-vector, it indicates the clipping
%		boundaries [min..max]
%
%
%EXAMPLES:
%       im(rand(5));            % show some image
%       im(rand(5),.4)          % clip the values to [0 .4]
%       im(rand(5), [.3 .4])    % clip the values to [.3 .4]
%
%OUT:
%	h - handle to the image (optional)

% $Id: im.m,v 1.15 2003/03/24 11:03:24 wexler Exp $
% Yoni, Fri Sep 28 15:03:33 2001

clip = length(varargin)>0;
if clip & length(varargin{1})==2
   clip_min=varargin{1}(1);
   clip_max=varargin{1}(2);
elseif clip & length(varargin{1})==1
   clip_min=0;
   clip_max=varargin{1};
else
   clip_min=0;
   clip_max=1;
end

if isstruct(x)
   if isfield(x,'img')
      x=x.img;
   else
      error('When given a struct, im() expects a field named ''img''');
   end
end

if ndims(x)>3                           % we have a movie?
   x=x(:,:,:,1);  % show first frame
end

if isa(x, 'double')
%   if any(isnan(x(:))) % This makes matlab display it much faster
%      m=min(x(:));
%      x(isnan(x))=m;
%   end
   if length(size(x)) == 3 & (min(x(:))<0 | max(x(:))>1)
      % 3-channel image, scale to 0..1

      if clip
         x(x<clip_min)=clip_min;
         x(x>clip_max)=clip_max;
      end

      m=min(x(:)); M=max(x(:));
      if(m<M)
	 x = (x-m)/(M-m);
      else
	 x = x-m;
      end
   end
   
end

if isa(x, 'logical')
    x=double(x);
end

%iptsetpref('ImshowBorder', 'tight')
hand=imagesc(x);
%cg(64);
colormap gray
axis image

if clip
   set(get(hand, 'parent'), 'ClimMode', 'manual', 'Clim', [clip_min clip_max]);
end

%set(gcbf, 'KeyPressFcn', 'zoom')

if nargout>0
   h=hand;
end

%pixval
%impixelinfo

