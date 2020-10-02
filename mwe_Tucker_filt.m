clear; close all; clc
rootDir = return_repository_root();
%addpath(genpath(fullfile(rootDir, '_toolboxes', 'xUnit', 'src')))
addpath(genpath(fullfile(rootDir, '_toolboxes', 'tensorlab_2016-03-28')))
addpath(genpath(fullfile(rootDir, '_toolboxes', 'Tensor_Toolbox_CRL_beta')))
addpath(genpath(fullfile(rootDir, '_utils')))


rng('default'); rng(7); %<-- Set seed for reproducibility

%% Setup and Parameters

% Filter dimensions and Rank
I = [4, 4, 4];
R = [3, 3, 3];
opts.max_it = 5;
N = length(I);

% System identification setup

L = prod(I); % number of channel taps
K = L; % sample size
sn = 12; % SNR in dB

%% Channel
h = (randn(L, 1) + 1i * randn(L, 1)) * sqrt(2);
%h = (1:L)';

%% Signal

% ZMCSCG signal with variance 1
x = (randn(K+L-1, 1) + 1i * randn(K+L-1, 1)) / sqrt(2);
pow = x' * x / (K+L-1); % should be rougly 1 (OK)

% Additive noise
sigma_n = sqrt(pow*10^(-sn / 10));
noi = sigma_n * (randn(K, 1) + 1i * randn(K, 1)) * sqrt(2);

% Observed signal
H = convmtx(h, K+L-1);
d = H(L:K+L-1, :) * x + noi; % skip L - 1 samples in the beginning...


%% SO-Statistics

% Prepare Toeplitz matrix, as in lecture
X = nan(L, K);
for k = 1:K
    X(:, k) = x(k+L-1:-1:k);
end


%% Wiener-Hopf Solution
if(1)
    tic
    R_hat = 1 / K * (X * X');
    p_hat = 1 / K * X * conj(d);
    
    w_opt = R_hat \ p_hat;
    toc
    err_est_WH = norm(h - conj(w_opt), 2)^2 ./ norm(conj(h), 2)^2
    
    if sigma_n == 0
        assertElementsAlmostEqual(h, conj(w_opt)) % OK
    end
end
%% Tucker-filter
%R = mlrank(reshape(X(:, 1), I), 1e-5)
%mlrankest(reshape(X, [I, K]));


Gt = [];
Gt = supereye(N, 3);
tic
w_opt_Tuck = Tucker_filt(X, d, I, R, opts, Gt);
toc

err_est_Tuck = norm(h - conj(w_opt_Tuck), 2)^2 ./ norm(conj(h), 2)^2



figure
plot(abs(h), '-r', 'Linewidth', 0.5); hold on, grid on;
plot(abs(w_opt_Tuck), '-.b', 'Linewidth', 1);
legend('Unknown Channel', 'Estimated Channel');
