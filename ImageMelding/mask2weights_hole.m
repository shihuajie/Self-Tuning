function weights = mask2weights_hole(mask, alpha_weight)
%ObjInfluenceDistance = 0.33*nnz(mask)^0.5; % the influence region of an object is proportional to it's size, assuming one object in mask
%ObjInfluenceDistance = 0.05*numel(mask)^0.5; % the influence region of an object is proportional to the image size
disttrans = bwdist(mask);
%disttrans(disttrans>ObjInfluenceDistance) = ObjInfluenceDistance;
weights = 5+200*alpha_weight.^(-disttrans-5*double(~mask));
end