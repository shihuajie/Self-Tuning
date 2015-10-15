function hfig = fim(varargin)

hfig = [];
if isempty(varargin{1}) || (~isa(varargin{1},'double') && ~isa(varargin{1},'single') && ~isa(varargin{1},'uint8') && ~isa(varargin{1},'int32') && ~isa(varargin{1},'logical'))
	disp('fim: Not and image or image format non-standard!') 
	return
end
%hfig = figure('KeyPressFcn', @printfig); 
if size(varargin{1},1)==1 && size(varargin{1},2)==1
	hfig = figure(varargin{1});
	set(hfig,'Color',[0.8 0.8 0.8], ...
		'Name','im', ...
		'KeyPressFcn', @printfig, ...
		'Pointer', 'crosshair', ...
		'BusyAction', 'queue', ...
		'DoubleBuffer', 'on', ...
		'Units','normalized');
else
	hfig = figure('Color',[0.8 0.8 0.8], ...
		'Name','im', ...
		'KeyPressFcn', @printfig, ...
		'Pointer', 'crosshair', ...
		'BusyAction', 'queue', ...
		'DoubleBuffer', 'on', ...
		'Units','normalized');
% 	, ...
%	'NumberTitle','on', ...
% 	'ButtonDownFcn', 'disp(''Click on images'')',...
% 	'CloseRequestFcn', 'fim close', ...
end

if size(varargin{1},1)==1 && size(varargin{1},2)==1
	ime(varargin{2},varargin{3:end});
else
	ime(varargin{1},varargin{2:end});
end

%pixval
impixelinfo



function [] = printfig(src, evnt)

hfig = gcf;
%switch get(hfig, 'CurrentCharacter')
switch evnt.Character
	case {'m'}
		pos=get(hfig, 'Position');
		dx=pos(3)/10; dy=pos(4)/10;
		set(hfig, 'Position', pos+[-dx -dy 2*dx 2*dy]);
	case {'b'}
		pos=get(hfig, 'Position');
		dx=pos(3)/10; dy=pos(4)/10;
		set(hfig, 'Position', pos-[-dx -dy 2*dx 2*dy]);
	case {'t','T','n','N'}
		truesize
	case {'j'}
		colormap jet;
		%colorbar;
	case {'h'}
		colormap hsv;
		%colorbar;
	case {'g'}
		colormap gray;
	case {'d'}
		imdistline;
		
end

