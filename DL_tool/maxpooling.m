function [res, maxpos] = maxpooling(array,scale)
% @ array: input data
% @ scale: pooling factor
% @ res: max_pooling value
% @ maxpos: position of max_pooling value
    num = floor( size( array(:,:,1) )/scale );
    for i = 1:size(array,3)
        for j = 1 : num(2)
            res(1,j,i) = max(array(1, (j-1)*scale + 1 : j * scale, i));
            maxpos(1,j,i) = find(array(1, (j-1)*scale + 1 : j * scale, i) == res(1,j,i),1,'first') + (j-1)*scale;
        end
    end
end