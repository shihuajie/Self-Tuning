%%
% smart_init - initialize a target image using random patches smartly
%
%% Input
%   trg_im      the target image
%   src_im      the source texture
%   G           the generators [g1; g2]
%   P           the patch size
%   alignment_channels
%               the channels to use for alignment purpose
%   method      the alignment method from
%                   - blocks: use neighboring blocks to find offsets
%                   - block_boundaries: use neighboring block boundaries
%                   - pixels: pixel-based offsetting
%   trys        the number of block attempts
%
%% Output
%   out_im      an image of the target size, with randomly
%               shuffled patches from the source texture
function [ out_im ] = smart_init( trg_im, src_im, G, options )

    if nargin < 4
        options = struct([]);
    end
       
    % workspace
    [trgSzX, trgSzY, nCh] = size(trg_im);
    out_im = zeros(trgSzX, trgSzY, nCh);
    mask = logical(out_im(:, :, 1)); % the areas which have been processed
    [srcSzX, srcSzY, ~] = size(src_im);

    % default parameters
    P = get_option(options, 'patch_size', 10);
    if nCh >= 7
        alignment_channels = get_option(options, 'alignment_channels', [1,7]);
    else
        alignment_channels = get_option(options, 'alignment_channels', 1);
    end
    method = get_option(options, 'smart_init_method', 'block_boundaries');
    trys = get_option(options, 'smart_init_trys', 5);
    
    if strcmp(method, 'block_boundaries')
        method = 'blocks';
        use_boundaries = 1;
    else
        use_boundaries = 0;
    end
    
    switch lower(method)
        case 'blocks'
          %% start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % johnsback 2015-08-03
            %%read source & target guidance channel
%             src_gc_file = 'src_gc.png';
%             trg_gc_file = 'trg_gc.png'
%             nnf_file = 'test_gc_nnf.bmp';
%             nnf_matfile = 'test_gc_nnf.mat';
%             global out_directory;
%             if ~isempty(out_directory)
%                 src_gc_file = strcat(out_directory, '\', src_gc_file);
%                 trg_gc_file = strcat(out_directory, '\', trg_gc_file);
%                 nnf_file = strcat(out_directory, '\', nnf_file);
%                 nnf_matfile = strcat(out_directory, '\', nnf_matfile);
%             end
%             src_gc = imread(src_gc_file);
%             trg_gc =  imread(trg_gc_file);
% 
%             % do patch match for initialization
% %             if(~exist(nnf_matfile, 'file'))
% %                 testnnf = [];
% %                 testnnf_msk = [];
% %                 [trg_NNF, dist] = nnmex(trg_gc, src_gc, testnnf, testnnf_msk, options);
% %                 save(nnf_matfile, 'trg_NNF', '-ascii');
% %             else
% %                 trg_NNF = load(nnf_matfile, '-ascii');               
% %             end
%                             
%             trg_gc_r = trg_gc(:, :, 1);
%             trg_gc_g = trg_gc(:, :, 2);
%             trg_gc_b = trg_gc(:, :, 3);
%             trg_gc_rg = trg_gc_r;% + trg_gc_g;
%             thresh = graythresh(trg_gc_rg);
%             trg_bin_msk = im2bw(trg_gc_rg, thresh);
%            %% imwrite(trg_bin_msk, strcat(out_directory, '\test_bin.png'));
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%             % Initialization by pixel            
%             srcPos_blue_n = 0; srcPos_yellow_n = 0; srcPos_pink_n = 0;
%             srcPos_blue = zeros(srcSzX*srcSzY, 2);
%             srcPos_yellow =zeros(srcSzX*srcSzY, 2);
%             srcPos_pink = zeros(srcSzX*srcSzY, 2);
%             for xx = 1 :  srcSzX
%                 for yy = 1 : srcSzY
%                     if(src_gc(xx, yy, 1) > 128) % red
%                         if(src_gc(xx, yy, 2) > 128) % yellow
%                             srcPos_yellow_n = srcPos_yellow_n + 1;
%                             srcPos_yellow(srcPos_yellow_n, :) = [xx, yy];
%                         else
%                             srcPos_pink_n = srcPos_pink_n + 1;
%                             srcPos_pink(srcPos_pink_n, :) = [xx, yy];
%                         end
%                     else
%                         srcPos_blue_n = srcPos_blue_n + 1;
%                         srcPos_blue(srcPos_blue_n, :) = [xx, yy];
%                     end
%                 end
%             end
%  
%             randPos = rand(trgSzX, trgSzY);
%             for xx = 1 :  trgSzX
%                 for yy = 1 : trgSzY
%                     randIdx = randPos(xx, yy);
%                      if(trg_gc(xx, yy, 1) > 128) % red
%                         if(trg_gc(xx, yy, 2) > 128) % yellow
%                             posIdx =  ceil(randIdx*srcPos_yellow_n);
%                             srcXX = srcPos_yellow(posIdx, 1);
%                             srcYY = srcPos_yellow(posIdx, 2);
%                         else
%                             posIdx = ceil(randIdx*srcPos_pink_n);
%                             srcXX = srcPos_pink(posIdx, 1);
%                             srcYY = srcPos_pink(posIdx, 2);
%                         end;
%                      else if(trg_gc(xx, yy, 2) > 128) %green
%                             srcXX = ceil(rand()*srcSzX);
%                             srcYY = ceil(rand()*srcSzY);
%                          else
%                              posIdx = ceil(randIdx*srcPos_blue_n);
%                              srcXX = srcPos_blue(posIdx, 1);
%                              srcYY = srcPos_blue(posIdx, 2);
%                          end
%                      end;
%                     out_im(xx, yy, :) = src_im(srcXX, srcYY, :);
%                 end
%             end
%          out_im(69:69+374, 7:506, :) = src_im(1:375, 1:500, :);


%% added by huajie 2015-9-15
%% start, initial by label to label
            srcPos_lab = zeros(srcSzX*srcSzY, 27, 2,'int32');
            num_lab = zeros(27,1,'int32');
            for xx = 1:srcSzX
                for yy = 1:srcSzY
                    label = round(src_im(xx,yy,8))+1;
                    index = num_lab(label) + 1;
                    srcPos_lab(index, label, 1) = xx;
                    srcPos_lab(index, label, 2) = yy;
                    num_lab(label) = index;
                end
            end
            randPos = rand(trgSzX, trgSzY);      
            for xx = 1:trgSzX
                for yy = 1:trgSzY
                    label = int32(round(trg_im(xx,yy,8))+1);
                    if(label == 1 || label == 27)
                        xPos = ceil(rand()*srcSzX);
                        yPos = ceil(rand()*srcSzY);
                    else
                        randIdx = ceil(randPos(xx, yy)*num_lab(label));
                        if(int32(randIdx) == 0)
                            continue;
                        end
                        xPos = srcPos_lab(randIdx, label, 1);
                        yPos = srcPos_lab(randIdx, label, 2);
                    end
                    out_im(xx, yy, 1:7) = src_im(xPos, yPos, 1:7);
                end
            end
            out_im(:, :, 8) = trg_im(:, :, 8);
%% end
          return;
% %% end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % use the lattice generators if any, or the patch size
            P_min = min([G; 0, 0], [], 1);
            P_max = max([G; 0, 0], [], 1);
            P_x = max(P_max(1) - P_min(1), 2 * P);
            P_y = max(P_max(2) - P_min(2), 2 * P);

            % split the block sizes if they're too big
            while trgSzX / P_x <= 3 || srcSzX / P_x <= 3
                P_x = ceil(P_x / 2);
            end
            while trgSzY / P_y <= 3 || srcSzY / P_y <= 3
                P_y = ceil(P_y / 2);
            end
            
            % boundary margin
            if use_boundaries
                % boundary margin
                bRX = get_option(options, 'smart_init_margin', 25); % ceil(P / 2));
                bRY = get_option(options, 'smart_init_margin', 25); % ceil(P / 2));
                % limited by the size of P_x and P_y
                bRX = min(bRX, P_x - 1);
                bRY = min(bRY, P_y - 1);
            end

            % grid of blocks
            trgDimX = ceil(trgSzX / P_x);
            trgDimY = ceil(trgSzY / P_y);
            srcDimX = floor(srcSzX / P_x);
            srcDimY = floor(srcSzY / P_y);

            % blocks positions
            trgCount = trgDimX * trgDimY;
            srcCount = srcDimX * srcDimY;
            totalCount = trgCount * trys;
            Iperm = [];
            while numel(Iperm) < totalCount
                Iperm = [Iperm randperm(srcCount)];
            end
            Iperm = Iperm(1:totalCount);
            % reshape
            [IpermX, IpermY] = ind2sub([srcDimX, srcDimY], Iperm);
            posX = reshape(IpermX, [trgDimX, trgDimY, trys]) - 1; % 0-based indexing
            posY = reshape(IpermY, [trgDimX, trgDimY, trys]) - 1;
            % randi(srcDimX - 1, trgDimX, trgDimY, trys);
            % randi(srcDimY - 1, trgDimX, trgDimY, trys);

            % build the puzzle, row by row, block by block
            for i = 1 : trgDimX
                % beginning and end of the block (X)
                dstXS = (i - 1) * P_x + 1;
                dstXE = min(i * P_x, trgSzX);
                lX = dstXE - dstXS; % effective size of the block (X)

                for j = 1 : trgDimY
                    fprintf('Block (%d,%d)\n', i, j);
                    % beginning and end of the block (Y)
                    dstYS = (j - 1) * P_y + 1;
                    dstYE = min(j * P_y, trgSzY);
                    lY = dstYE - dstYS; % effective size of the block (Y)
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Initialization by block
                    count = 0;
                    for jj = dstYS:dstYE
                        for ii = dstXS:dstXE
                            if(trg_bin_msk(ii, jj) > 0)
                                count = count +1;
                            end
                        end
                    end
                    
                    if(count > 5) % use the patch match result
                        % vx = trg_NNF(dstXS, dstYS, 1);
                        % vy = trg_NNF(dstXS, dstYS, 2);
                        
                        bestfit_block = 9999999999;
                        bestfit_block_xx = 0;
                        bestfit_block_yy = 0;
                        for xx = 1:(srcSzX - lX)
                            for yy = 1:(srcSzY - lY)
                                trg_block = trg_gc(dstXS : dstXE, dstYS : dstYE, :);
                                src_block = src_gc(xx : xx+lX, yy : yy+lY, :);
                                dd_block = trg_block - src_block;
                                dd = sum(sum(sum(dd_block.*dd_block)));
                                if(dd < bestfit_block)
                                    bestfit_block = dd;
                                    bestfit_block_xx = xx;
                                    bestfit_block_yy = yy;
                                end
                            end
                        end
                        
                        out_im(dstXS : dstXE, dstYS : dstYE, :) = src_im(bestfit_block_xx : bestfit_block_xx+lX, bestfit_block_yy : bestfit_block_yy+lY,:);
                        mask(dstXS : dstXE, dstYS : dstYE) = 1; % fill the mask
                        continue;
                    end
                                   
                   continue;

                   % Initialization by block
%% %modified by huajie start
%                     count_yellow = 0;
%                     count_pink = 0;
%                     count_blue = 0;
%                     for ii = dstYS:dstYE
%                         for jj = dstXS:dstXE
% 	  						if trg_gc(ii, jj, 1) > 128 % red
% 		                        if trg_gc(ii, jj, 2) > 128 % yellow
% 		                            count_yellow = count_yellow + 1;
% 		                        else
% 		                            count_pink = count_pink + 1;
% 		                        end
% 		                    else
% 		                        count_blue = count_blue + 1;
% 		                    end
%                         end
%                     end
% 
% %                     if count_blue > count_pink
% %                     	if count_blue > count_yellow
% %                     		maxColor = 'blue';
% %                     	else 
% %                     		maxColor = 'yellow';
% %                     	end
% %                     else
%                     	if count_pink > count_yellow
%                     		maxColor = 'pink';
%                     	else
%                     		maxColor = 'yellow';
%                     	end
% %                     end
%                     bestfit_block = 9999999999;
%                     bestfit_block_xx = 0;
%                     bestfit_block_yy = 0;
% 
%                     trg_block = trg_gc(dstXS : dstXE, dstYS : dstYE, :);
%                     for xx = 1:(srcSzX - lX)
%                         for yy = 1:(srcSzY - lY)
%                             src_block = src_gc(xx : xx+lX, yy : yy+lY, :);
%                             dd_block = trg_block - src_block;
%                         	switch maxColor
%                         		case 'yellow'
%                          			dd = sum(sum(dd_block(:,:,2).*dd_block(:,:,2)));
%                         		case 'pink'
%                         			dd = sum(sum(dd_block(:,:,1).*dd_block(:,:,1)));
%                         		%case 'blue'
%                         		%	dd = sum(sum(dd_block(:,:,3).*dd_block(:,:,3)));
%                         	end
% 
%                             if(dd < bestfit_block)
%                                 bestfit_block = dd;
%                                 bestfit_block_xx = xx;
%                                 bestfit_block_yy = yy;
%                             end
%                         end
%                     end
%                     
%                     out_im(dstXS : dstXE, dstYS : dstYE, :) = src_im(bestfit_block_xx : bestfit_block_xx+lX, bestfit_block_yy : bestfit_block_yy+lY,:);
%                     mask(dstXS : dstXE, dstYS : dstYE) = 1; % fill the mask
%                                 
%                    continue;
%% % end ---huajie
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
                    % try a few blocks
                    bestK = -inf;
                    for t = 1:trys
                    
                        % beginning of the block in the source
                        srcXS0 = 1 + posX(i, j, t) * P_x; % this is only correct without noise
                        srcYS0 = 1 + posY(i, j, t) * P_y;
                        % - working with noise, we use correlation to find the best offset
                        if i * j > 1
                            % basic patch
                            zeroPatch = zeros(P_x * 3, P_y * 3);
                            num_align_ch = numel(alignment_channels);
                            
                            if isfield(options, 'debug') && options.debug == -43
                                srcXS = srcXS0;
                                srcYS = srcYS0;
                                break;
                            end

                            % the exemplar patch to match with xcorr
                            srcRgX = max(1, srcXS0 - P_x):min(srcSzX, srcXS0 + 2 * P_x - 1);
                            srcRgY = max(1, srcYS0 - P_y):min(srcSzY, srcYS0 + 2 * P_y - 1);
                            srcPatch = repmat(zeroPatch, [1, 1, num_align_ch]);
                            srcMask = logical(zeroPatch);
                            srcShX = P_x - srcXS0 + 1;
                            srcShY = P_y - srcYS0 + 1;
                            srcPatch(srcRgX + srcShX, srcRgY + srcShY, :) = src_im(srcRgX, srcRgY, alignment_channels);
                            srcMask(srcRgX + srcShX, srcRgY + srcShY) = 1;

                            % the target patch
                            trgRgX = max(1, dstXS - P_x):min(trgSzX, dstXS + P_x - 1);
                            trgRgY = max(1, dstYS - P_y):min(trgSzY, dstYS + 2 * P_y - 1);
                            trgPatch = repmat(zeroPatch, [1, 1, num_align_ch]);
                            trgMask = logical(zeroPatch);
                            trgShX = P_x - dstXS + 1;
                            trgShY = P_y - dstYS + 1;
                            trgPatch(trgRgX + trgShX, trgRgY + trgShY, :) = out_im(trgRgX, trgRgY, alignment_channels);
                            trgMask(trgRgX + trgShX, trgRgY + trgShY) = mask(trgRgX, trgRgY);

                            % boundary mode
                            if use_boundaries
                                boundaryMask = logical(zeroPatch);
                                boundaryMask((P_x - bRX):(P_x + bRX), (P_y - bRY):(2 * P_y + bRY)) = 1;
                                boundaryMask((2 * P_x - bRX):(2 * P_x + bRX), (P_y - bRY):(2 * P_y + bRY)) = 1;
                                boundaryMask((P_x - bRX):( 2 * P_x + bRX), (P_y - bRY):(P_y + bRY)) = 1;
                                boundaryMask((P_x - bRX):( 2 * P_x + bRX), (2 * P_y - bRY):(2 * P_y + bRY)) = 1;
                                trgMask = trgMask & boundaryMask; % only correlation with boundary regions
                            end

                            % keep only the parts that we can use during correlation
                            minShX = srcXS0 - 1;
                            minShY = srcYS0 - 1;
                            maxShX = srcSzX - lX - srcXS0; % - 1;
                            maxShY = srcSzY - lY - srcYS0; % - 1;

                            % use correlation to find best translation
                            % subplot(1, 2, 1); imshow(srcPatch); subplot(1, 2, 2); imshow(trgPatch); drawnow;
                            [tX, tY, K] = find_alignment(trgPatch, srcPatch, trgMask, srcMask, ...
                                @(tx,ty) and(...
                                    and(-tx <= maxShX, -ty <= maxShY), ...
                                    and(-tx >= -minShX, -ty >= -minShY) ...
                                ), ... % validity
                                1 ... % fast method
                            );
                            fprintf('tX=%g, tY=%g, K=%g\n', tX, tY, K);

                            % best shift
                            if K > bestK
                                bestK = K;
                                srcXS = srcXS0 - tX;
                                srcYS = srcYS0 - tY; % not +tY!!!
                            end
                        else
                            srcXS = srcXS0;
                            srcYS = srcYS0;
                            break; % no need to try more at the beginning 
                        end
                    end

                    % prevent death
                    if srcXS <= 0
                        fprintf('srcXS--, srcXS=%d\n', srcXS);
                        srcXS = 1;
                    elseif srcXS + lX > srcSzX
                        fprintf('srcXS++, srcXS=%d\n', srcXS);
                        srcXS = srcSzX - lX;
                    end
                    if srcYS <= 0
                        fprintf('srcYS--, srcYS=%d\n', srcYS);
                        srcYS = 1;
                    elseif srcYS + lY > srcSzY
                        fprintf('srcYS++, srcYS=%d\n', srcYS);
                        srcYS = srcSzY - lY;
                    end

                    % use that patch
                    out_im(dstXS : dstXE, dstYS : dstYE, :) = src_im(srcXS : srcXS+lX, srcYS : srcYS+lY,:);
                    mask(dstXS : dstXE, dstYS : dstYE) = 1; % fill the mask
                    if isfield(options, 'debug') && options.debug == -41
                        keyboard
                    end
                    % subplot(1, 2, 1); imshow(out_im(:, :, 1) / 100);
                    % subplot(1, 2, 2); imshow(out_im(:, :, 7));
                    % drawnow;
                end
                
            end
            if isfield(options, 'debug') && options.debug == -42
               keyboard 
            end
        
        case 'pixels'
            g1 = G(1, 1:2);
            g2 = G(2, 1:2);
            
            % pixels locations
            [X0, Y0] = ndgrid(1:trgSzX, 1:trgSzY);
            % wrap into coordinates within the source
            X = mod(X0, srcSzX) + 1;
            Y = mod(Y0, srcSzY) + 1;
            % use a * g1 + b * g2
            maxA = maxMultiple([srcSzX, srcSzY], g1)
            maxB = maxMultiple([srcSzX, srcSzY], g2)
            
            A = randi([0, maxA], trgSzX, trgSzY);
            B = randi([0, maxB], trgSzX, trgSzY);
            
            X = mod(X + A * g1(1) + B * g2(1) + 10 * srcSzX, srcSzX) + 1;
            Y = mod(Y + A * g1(2) + B * g2(2) + 10 * srcSzY, srcSzY) + 1;
            % TODO do not use mod as we don't know if the repetition is
            % aligned with the image size!!!
            I = sub2ind([srcSzX, srcSzY], X(:), Y(:));
            
            % work on each channel
            for c = 1:nCh
                src_ch = src_im(:, :, c);
                out_im(:, :, c) = reshape(src_ch(I), trgSzX, trgSzY);
            end
    end

end

function value = get_option(options, name, default_value)
    if isstruct(options) && isfield(options, name)
        value = options.(name);
    else
        value = default_value;
    end
end

function m = maxMultiple(sz, g)
    m = 1;
    for i = 1:numel(sz)
        v = ceil(sz(i) / abs(g(i)));
        if ~isinf(v) && ~isnan(v)
           m = max(m, v);
        end
    end
end