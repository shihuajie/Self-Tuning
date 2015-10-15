function [ h ] = multimshow( grid, varargin )
%MULTIMSHOW Show multiple images in a big subplot figure
    nX = grid(1);
    if numel(grid) > 1
        nY = grid(2);
        start = 1;
    else
        nY = varargin{1};
        start = 2;
    end
    
    method=@imshow;
    
    idx = 1;
    for i = start:numel(varargin)
        data = varargin{i};
        if isnumeric(data) || islogical(data)
            subplot(nX, nY, idx);
            method(data);
            idx = idx+1;
        elseif isa(data, 'function_handle')
            method=data;
        end
    end
end

