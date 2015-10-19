function [output_NNF] = rescale_NNF( cur_NNF, cur_trg_sz, cur_tex_sz, prev_tex_sz, options )

% the data we're working with
cur_trg_szX = cur_trg_sz(1);
cur_trg_szY = cur_trg_sz(2);
cur_tex_szX = cur_tex_sz(1);
cur_tex_szY = cur_tex_sz(2);
prev_tex_szX = prev_tex_sz(1);
prev_tex_szY = prev_tex_sz(2);

%% Basic data resampling
szX = size(cur_NNF, 1);
szY = size(cur_NNF, 2);
szC = size(cur_NNF, 3);
P = options.patch_size;
S = P - 1;
% resize the nnf using the correct image size (not the shorter nnf size)
padded_NNF = zeros(szX + S, szY + S, szC, class(cur_NNF));
padded_NNF(1:szX, 1:szY, :) = cur_NNF;
padded_NNF(1:szX, (szY+1):end, :) = repmat(cur_NNF(1:szX, szY, :), 1, S);
padded_NNF((szX+1):end, 1:szY, :) = repmat(cur_NNF(szX, 1:szY, :), S, 1);
padded_NNF((szX+1):end, (szY+1):end, :) = repmat(cur_NNF(szX, szY, :), S, S);
% upsampling with the same ratio as the target image
output_NNF = imresize(padded_NNF, cur_trg_sz, 'bilinear'); %% modified by huajie 2015-9-25 bilinear to lanczos3
output_NNF = output_NNF(1:(cur_trg_szX - S), 1:(cur_trg_szY - S), :);

%% Basic patch field rescaling
% rescale the Y+X channels
rtX = cur_tex_szX / prev_tex_szX;
rtY = cur_tex_szY / prev_tex_szY;
% ... and boundary checking
sX = rtX * output_NNF(:,:,1);
sX(sX < 0) = 0;
sX(sX > cur_tex_szX - P) = cur_tex_szX - P; % /!\ nnf is 0-based!
sY = rtY * output_NNF(:,:,2);
sY(sY < 0) = 0;
sY(sY > cur_tex_szY - P) = cur_tex_szY - P;
% output
output_NNF(:, :, 1) = sX;
output_NNF(:, :, 2) = sY;

%% Affine patch field rescaling
if size(cur_NNF, 3) > 2
    % rotation bound checking
    rt = output_NNF(:,:,3);
%   rt(rt < options.min_angle) = options.min_angle; %% modified by huajie 2015-10-15
%   rt(rt > options.max_angle) = options.max_angle; %% index no need to check
    % scaling bound checking
    scX = output_NNF(:,:,4);
    scX(abs(scX) < options.min_scale) = options.min_scale * sign(scX(abs(scX) < options.min_scale));
    scX(abs(scX) > options.max_scale) = options.max_scale * sign(scX(abs(scX) > options.max_scale));
    scY = output_NNF(:,:,5);
    scY(abs(scY) < options.min_scale) = options.min_scale * sign(scY(abs(scY) < options.min_scale));
    scY(abs(scY) > options.max_scale) = options.max_scale * sign(scY(abs(scY) > options.max_scale));
    % output
    output_NNF(:, :, 3) = rt;
    output_NNF(:, :, 4) = scX;
    output_NNF(:, :, 5) = scY;
end
