function X = iunfolding(Xn, n, I, order)
% Reconstructs a tensor out of its n-mode unfolding.
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
% |     Last modifications: 03.07.2019
% |         - minor edits on notation
% |     13.06.2019 (MG)
% |         -elaborated on description of column orderings
% |     Date authored: 04.23.2007
% |----------------------------------------------------------------
%
% X = iunfolding(Xn, n, I, order)
%
% reproduces the origin tensor X out of its n-mode unfolding
% Matrix Xn, produced by the function unfolding (type 'help unfolding'
% to get additional informations). Therefore, the dimensions of the
% origin tensor have to be given by the vector sizes. The optional
% parameter order is the ordering used by the unfolding command. If order
% is not given, the function assumes order = 3 (Lathauwer unfolding).
% Otherwise the following ordering of the columns is assumed:
%
%   order = 1: forward column ordering (MATLAB, column-major)
%       indices of n - mode vectors go faster with increasing index
%   order = 2: reverse column ordering (Tensorly, row-major)
%       indices of n - mode vectors go slower with increasing index
%   order = 3: reverse cyclical, bc (de Lathauwer)
%       indices go slower with I_n+1, ... I_N, I_1, ... I_n-1
%   order = 4: forward cyclical, fc (flipped de Lathauwer)
%       indices go slower with I_n-1, ... I_1, I_N, ... I_n+1
%
% Inputs: Tn    - matrix with the n - mode vectors of a tensor T
%         n     - dimension
%         sizes - vector containg the size of T
%         order - defines the ordering of the n - mode vectors (optional)
%
% Output: T     - reproduced tensor

% get dimension
N = length(I);

% make singletons at the end of T possible
if n > N
    I = [I, ones(1, n-N)];
    N = n;
end

% Set standard Lathauwer unfolding
if nargin == 3
    order = 3;
end

% get permutation vector
switch order
    case 1
        permute_vec = [n, 1:(n - 1), (n + 1):N]; % indices go faster with increasing index
        [temp, ipermute_vec] = sort(permute_vec);
    case 2
        permute_vec = [n, fliplr([1:(n - 1), (n + 1):N])]; % indices go slower with increasing index
        [temp, ipermute_vec] = sort(permute_vec);
    case 3
        permute_vec = [n, fliplr(1:(n - 1)), fliplr((n + 1):N)]; % Lathauwer: indices go slower with I_n+1, ... I_N, I_1, ... I_n-1
        [temp, ipermute_vec] = sort(permute_vec);
    case 4
        permute_vec = [n, fliplr([fliplr(1:(n - 1)), fliplr((n + 1):N)])]; % flipped Lathauwer
        [temp, ipermute_vec] = sort(permute_vec);
    otherwise
        disp('Error: unknown ordering for n--mode vectors');
        return
end

% get origin tensor
X = permute(reshape(Xn, I(permute_vec)), ipermute_vec);
