function x = vec(X)
% VEC   Computes the vectorized representation of tensor X in matlab
% columns ordering.
%
% |----------------------------------------------------------------
% | (C) 2018 TU Ilmenau, Communications Research Laboratory
% |
% |     Mikus Grasis
% |
% |     Advisor:
% |        Univ. Prof. Dr.-Ing. Martin Haardt
% |
% |     Last modifications: 28/06/2018
% |----------------------------------------------------------------
%
% Syntax:
%    x = vec(X)
%
% Input:
%    X - input tensor
%
% Output:
%    x - vectorized version

x = X(:);

end
