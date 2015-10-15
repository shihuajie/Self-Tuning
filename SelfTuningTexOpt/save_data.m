function [] = save_data( data, file_name, file_ext, output_dir )
%SAVE_DATA Save a number or string into a file
%
% SYNOPSIS
%   save_data(data, file_name [, file_ext, output_dir]) - save data into $file_name
%   save_data(file_name) - create an empty file named $file_name
%
% INPUT
%   - data        the string or number to save
%   - file_name   the filename to use
%   - output_dir  the output directory (global out_directory by default)
%   - file_ext    the file extension (.dat or nothing if in file_name)
%

    if isnumeric(data)
        str = num2str(data);
    else
        str = data;
    end
    
    if nargin == 1
        file_name = data;
        str = '';
        file_ext = '';
    elseif nargin < 3
        if strfind(file_name, '.')
            file_ext = '';
        else
            file_ext = '.dat';
        end
    end
    global out_directory;
    if nargin < 4
        if isempty(out_directory)
            output_dir = '.';
        else
            output_dir = out_directory;
        end
    end
    
    try
        fid = fopen(fullfile(output_dir, [file_name file_ext]), 'w');
        if fid >= 0
            fprintf(fid, '%s', str);
            fclose(fid);
        else
            fprintf('Error saving data to %s!\n', file_name);
        end
        
    catch ex
        getMessage(ex)
    end
end

