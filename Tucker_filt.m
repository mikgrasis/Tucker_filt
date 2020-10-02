function [w_opt, info] = Tucker_filt(X, d, I, R, opts, Gt)
% Tensor MMSE filter based on Tucker model.
%
% |----------------------------------------------------------------
% | (C) 2020 TU Ilmenau, Communications Research Laboratory
% |
% |
% |  ______   ______     __   __     ______     ______     ______        ______   ______     ______     __         ______
% | /\__  _\ /\  ___\   /\ "-.\ \   /\  ___\   /\  __ \   /\  == \      /\__  _\ /\  __ \   /\  __ \   /\ \       /\  ___\
% | \/_/\ \/ \ \  __\   \ \ \-.  \  \ \___  \  \ \ \/\ \  \ \  __<      \/_/\ \/ \ \ \/\ \  \ \ \/\ \  \ \ \____  \ \___  \
% |    \ \_\  \ \_____\  \ \_\\"\_\  \/\_____\  \ \_____\  \ \_\ \_\       \ \_\  \ \_____\  \ \_____\  \ \_____\  \/\_____\
% |     \/_/   \/_____/   \/_/ \/_/   \/_____/   \/_____/   \/_/ /_/        \/_/   \/_____/   \/_____/   \/_____/   \/_____/
% |
% |
% |     Mikus Grasis
% |
% |     Advisors:
% |         Univ. Prof. Dr.-Ing. Martin Haardt
% |         Prof. Andre Lima Ferrer de Almeida
% |
% |     Date authored: 29.09.2020
% |----------------------------------------------------------------
%
%   [w_opt_tucker] = Tucker_filt(X, d, I, R)
%
% Inputs:   X       - signal matrix of size prod(I) x K (Toeplitz)
%           d       - desired response of size K x 1
%           I       - tensor filtering dimensions
%           R       - multilinear ranks
%           opts.max_it - maximum number of iterations
%           opts.tol    - tolerance for convergence (change of norm)
%
% Output:   w_opt   - Tucker-MMSE filter coefficients of size prod(I) x 1
if nargin < 6
    Gt = [];
end
if nargin < 5
    opts = [];
end
opts.max_it = setparam(opts, 'max_it', 5);
opts.tol = setparam(opts, 'tol', 1e-10);
assert(isequal(size(X, 1), prod(I)), 'number of filter coefficients must match size of input matrix')

K = size(X, 2); % number of data points
N = length(I);  % tensor order

% Sample statistics
R_xx = 1 / K * (X * X');
p_xd = 1 / K * X * conj(d);

% Tensorize input data
Xt = cell(1, K);
for k = 1:K
    Xt{k} = reshape(X(:, k), I);
end

% Initialize filter coefficients and core tensor
W = cell(1, N);
if isreal(X)
    for n = 1:N
        W{n} = randn(I(n), R(n));
    end
    if isempty(Gt)
        Gt = randn(R);
    end
else
    for n = 1:N
        W{n} = (randn(I(n), R(n)) + 1i * randn(I(n), R(n))) / sqrt(2);
    end
    if isempty(Gt)
        Gt = randn(R) + 1i * randn(R);
    end
end

% Alternating minimization
num_it = 0;
not_converged = true;
while not_converged && le(num_it, opts.max_it)
    num_it = num_it + 1;
    for n = 1:N
        %%% Estimate SO-statistics for n-th factor
        X_hat_n = nan(I(n)*R(n), K);
        W_minus_n_unf = unfolding(Gt, n, 1) * kron_Nd(W, n, 1).';
        for k = 1:K
            X_hat_n(:, k) = vec(unfolding(Xt{k}, n, 1) * W_minus_n_unf');
        end

        R_hat_n = 1 / K * (X_hat_n * X_hat_n');
        p_hat_n = 1 / K * X_hat_n * conj(d);
        
        %%% Update coefficients of n-th mode
        W{n} = reshape(R_hat_n \ p_hat_n, [I(n), R(n)]);
    end
    %%% Estimate SO-statistics for core tensor
    W_kron = kron_Nd(W, [], 1);
    R_hat = W_kron' * R_xx * W_kron;
    p_hat = W_kron' * p_xd;
    
    %%% Update coefficients of core tensor
    Gt = reshape(R_hat \ p_hat, R);
    
    
%     if(0)
%         % Check for convergence
%         Wt = cpdgen(W);
%         delta = norm(Wt(:) - Wt_old(:));
%         
%         Wt_old = Wt;
%         num_it = num_it + 1;
%         
%         if le(delta, opts.tol)
%             not_converged = false;
%         end
%         
%         % Show some output
%         if opt.verbose
%             fprintf('[TMMSE] it %d/%d, err = %d\n', num_it, opts.num_it, e);
%         end
%     end
    
end
info.num_it = num_it;
w_opt = kron_Nd(W, [], 1) * vec(Gt);
end
