function [ data ] = load_mat( file_name, var_name )
%LOAD_MAT Load the first variable found in a mat file
%
% INPUT
%   - file_name   the mat filename
%   - var_name    the variable name to expect (or get an error if not found)
%
% OUTPUT
%   - data        the first variable found in the mat file
%

    s = load(file_name);
    if nargin == 2
        if var_name == 1
            data = s;
        elseif isfield(s, var_name)
            data = s.(var_name);
        else
            error('Variable "%s" not found in %s', var_name, file_name);
        end
    else
        f = fieldnames(s);
        F = numel(f);
        if F < 1
            error('Empty mat file %s', file_name);
        else
            data = s.(f{1});
        end
    end
end

