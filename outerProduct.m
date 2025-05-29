function C = outerProduct(A, B)  % version 5
% OUTERPRODUCT computes the outer product of two matrices
% downloaded from Jan's solution in matlab central:
% https://www.mathworks.com/matlabcentral/answers/445798-outer-product-of-two-rectangular-matrices

C = reshape(A(:) * B(:).', [size(A), size(B)]);

end