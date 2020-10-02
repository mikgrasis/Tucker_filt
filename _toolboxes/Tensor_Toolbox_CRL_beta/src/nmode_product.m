function X_out = nmode_product(X_in, A, n)
% Computes the n-mode product of a tensor and a matrix.
%
% |----------------------------------------------------------------
% | (C) 2006 TU Ilmenau, Communications Research Laboratory
% |
% |     Martin Weis
% |
% |     Advisors:
% |        Dipl.-Ing. Giovanni Del Galdo
% |        Univ. Prof. Dr.-Ing. Martin Haardt
% |
% |     Last modifications: 03.07.2019 (MG)
% |         - minor edits on notation
% |     Date authored: 06.20.2006
% |----------------------------------------------------------------
%
% X_out = nmode_product(X_in, A, n)
%
% Computes the n-mode product of the tensor X_in and the matrix A.
% This means that all n-mode vectors of X are multiplied from the left
% with A.
%
% Inputs: X_in  - tensor
%         A     - matrix
%         n     - dimension
%
% Output: X_out - X_out = (X_in  x_n  A)

% Compute new dimensions of X_out
new_size = size(X_in);
new_size(n) = size(A, 1);

% Compute n-mode product
X_out = iunfolding(A*unfolding(X_in, n), n, new_size);
