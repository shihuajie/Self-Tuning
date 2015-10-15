function [ lines ] = read_lines( file, list )
%READ_LINES - read ascii lines from a file
%
% INPUT
%   - file    the file to read the lines from
%   - list    the filter index to use on the final lines
%
% OUTPUT
%   - lines   the filtered lines
    num_lines = 0;
    try
        fid = fopen(file);
        if fid
            tline = fgetl(fid);
            while ischar(tline)
                num_lines = num_lines + 1;
                lines{num_lines, 1} = tline; % remove end
                tline = fgetl(fid);
            end
            
            % filter
            if nargin == 2
                if numel(list) > 1
                    lines = lines(list);
                else
                    lines = lines{list};
                end
            end
        end
    catch ex
        getReport(ex)
    end
    if fid
        fclose(fid);
    end
end

