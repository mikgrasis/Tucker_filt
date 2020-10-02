function [w_opt, info] = CP_filt(X, d, I, R, opts)
% Tensor MMSE filter based on CP-model.
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
% |     Date authored: 15.07.2020
% |----------------------------------------------------------------
%
%   [w_opt_cp] = CP_filt(X, d, I, R)
%
% Inputs:   X       - signal matrix of size prod(L) x K (Toeplitz)
%           d       - desired response of size K x 1
%           I       - tensor filtering dimensions
%           R       - CP model rank
%           opts.max_it - maximum number of iterations
%           opts.tol    - tolerance for convergence (change of norm)
%
% Output:   w_opt   - tensor-MMSE filter coefficients of size prod(L) x 1
if nargin < 5
    opts = [];
end
opts.max_it = setparam(opts, 'max_it', 50);
opts.tol = setparam(opts, 'tol', 1e-10);
assert(isequal(size(X, 1), prod(I)), 'number of filter coefficients must match size of input matrix')

K = size(X, 2); % number of data points
N = length(I);  % tensor order

% Tensorize input data
Xt = reshape(X, [I, K]);

% Initialize filter coefficients
W = cell(1, N);
if isreal(X)
    for n = 1:N
        W{n} = randn(I(n), R);
    end
else
    for n = 1:N
        W{n} = (randn(I(n), R) + 1i * randn(I(n), R)) / sqrt(2);
    end
end

% Alternating minimization
num_it = 0;
not_converged = true;
while not_converged && le(num_it, opts.max_it)
    num_it = num_it + 1;
    for n = 1:N
        %%% Estimate SO-statistics
        X_bar_n = nan(I(n)*R, K);
        all_but_n = 1:N;
        all_but_n(n) = [];
        for r = 1:R
            X_bar_n_r = Xt;
            for ell = all_but_n
                X_bar_n_r = nmode_product(X_bar_n_r, W{ell}(:, r)', ell);
            end
            idx_start = 1 + I(n) * (r - 1);
            idx_end = idx_start + I(n) - 1;
            X_bar_n(idx_start:idx_end, :) = squeeze(X_bar_n_r);
        end
        R_hat_n = 1 / K * (X_bar_n * X_bar_n');
        p_hat_n = 1 / K * X_bar_n * conj(d);
        
        %%% Update coefficients of n-th mode
        W{n} = reshape(R_hat_n \ p_hat_n, [I(n), R]);
    end
    
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
w_opt = vec(cpdgen(W));
end
