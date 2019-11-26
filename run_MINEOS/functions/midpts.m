function [ midpt_matrix ] = midpts( matrix )
% Find the midpoints along the rows of a matrix

nrow = size(matrix,1);
ncol = size(matrix,2)-1;
midpt_matrix = zeros(nrow,ncol);
for irow = 1:nrow
    midpt_matrix(irow,:) = [matrix(irow,1:end-1)+matrix(irow,2:end)]/2;
end

