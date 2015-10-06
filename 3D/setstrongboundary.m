function [ norm_matrix ] = setstrongboundary( matrix, boundary, norm )
%Will replace points on the boundary of the matrix with the gradient passed
%in.

%TODO: validate input

norm_matrix = zeros(size(matrix));
matrix = matrix .* ~boundary;

%We can eliminate this if we ensure norm only exists along the boundary. 
matrix = matrix + norm .* boundary;

norm_matrix = matrix;

end

