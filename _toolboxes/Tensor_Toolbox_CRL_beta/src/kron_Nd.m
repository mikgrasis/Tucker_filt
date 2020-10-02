function B = kron_Nd(A, skip_dims, reverse)
% KRON_ND   N-D Kronecker product
%
% |----------------------------------------------------------------
% | (C) 2019 TU Ilmenau, Communications Research Laboratory
% |
% |     Mikus Grasis
% |
% |     Advisors:
% |         Univ. Prof. Dr.-Ing. Martin Haardt
% |         Prof. Andre Lima Ferrer de Almeida
% |
% |     Date authored: 22.09.2019
% |-----------------------------------------------------------------
%
% B = kron_Nd(A)
%
% computes the Kronecker product of all matrices given in the length-N
% cell array A, viz.,
%           B = KRON(KRON(KRON(..., A{N}), A{N-1}), ..., A{2}), A{1}).
%
% Inputs: A         - cell array of N matrices
%         skip_dims - vector with dimensions that should be skipped
%         reverse   - switch to reverse order of Kronecker product
%                       default: false
%
% Output: B         - Kronecker product of matrices
if nargin < 3, reverse = false; end
if nargin < 2, skip_dims = []; end

N = length(A);

set_N = 1:N;
set_N(skip_dims) = [];  % remove dimensions to be skipped
if reverse
    set_N = fliplr(set_N);
end

B = A{set_N(1)};
for curr_n = 2:length(set_N)
    B = kron(B, A{set_N(curr_n)});
end
end
