function [output_NNF] = upscale_NNmap(cur_NN, ...
  cur_trg_szX, cur_trg_szY, ...
  prev_src_szX, prev_src_szY, ...
  cur_src_szX, cur_src_szY, options)
output_NNF = imresize(cur_NN, [cur_trg_szX, cur_trg_szY], 'nearest');
rtX = cur_src_szX / prev_src_szX;
rtY = cur_src_szY / prev_src_szY;
sX = rtX * output_NNF(:,:,1);
sX(sX < 0) = 0;
sX(sX > cur_src_szX - options.patch_size + 1) = cur_src_szX - options.patch_size + 1;
sY = rtY * output_NNF(:,:,2);
sY(sY < 0) = 0;
sY(sY > cur_src_szY - options.patch_size + 1) = cur_src_szY - options.patch_size + 1;
rt = output_NNF(:,:,3);
rt(rt < options.min_angle) = options.min_angle;
rt(rt > options.max_angle) = options.max_angle;
sc = output_NNF(:,:,4);
sc(abs(sc) < options.min_scale) = options.min_scale;
sc(abs(sc) > options.max_scale) = options.max_scale;
output_NNF(:,:,1) = sX;
output_NNF(:,:,2) = sY;
output_NNF(:,:,3) = rt;
output_NNF(:,:,4) = sc;

