function results = fit_SEIR_model_OOS_C(dataTable, start_date, end_date, N, weights, normalize, num_trials, train_days)
% fit_SEIR_model_OOS_C - Out-of-sample variant of fit_SEIR_model_C.
%
% STANDALONE file - does not modify fit_SEIR_model_C. Same train/test
% split logic as fit_SIR_model_OOS_C (see its comments for full
% explanation): fitting only on the first train_days days, AICc on that
% portion, RMSE_test on the remaining days of the window (never seen by
% fitting), via simulation of the already-fitted model without re-fitting.

MIN_TEST = 5;

% === Data preprocessing (identical to fit_SEIR_model_C) ===
dates     = datetime(dataTable.Date);
confirmed = dataTable.Confirmed_Cases_corrected;
cumulative_cases = dataTable.Total_Cases;

window   = 10;
infected = movsum(confirmed, [window-1 0]);
infected = infected(:);
infected(infected < 0) = 0;

removed = cumulative_cases - infected;
removed = removed(:);
removed(removed < 0) = 0;

idx          = (dates >= start_date & dates <= end_date);
dates_fit    = dates(idx);
infected_fit = smoothdata(infected(idx), 'movmean', 7);
removed_fit  = smoothdata(removed(idx),  'movmean', 7);

T = length(dates_fit);

if T < train_days + MIN_TEST
    warning('fit_SEIR_model_OOS:tooShort', ...
        'Window %s-%s: T=%d, insufficient for train_days=%d (needs >=%d test points). Skipped.', ...
        datestr(start_date,'yyyy-mm-dd'), datestr(end_date,'yyyy-mm-dd'), T, train_days, MIN_TEST);
    results = [];
    return;
end

T_train = train_days;
T_test  = T - T_train;

infected_train = infected_fit(1:T_train);
removed_train  = removed_fit(1:T_train);

% === Initial conditions: I0, R0_init from the data; E0=I0 (assumption, same as fit_SEIR_model_C) ===
I0      = max(infected_fit(1), 1);
R0_init = removed_fit(1);
E0      = I0;
S0      = N - I0 - R0_init - E0;

assert(S0 > 0, 'S0 <= 0: check N or initial conditions.');

% === Objective function: ONLY over the training portion ===
wI = weights(1); wR = weights(2);
objective = @(params) compute_error(params, S0, E0, I0, R0_init, N, T_train, ...
                                    infected_train, removed_train, wI, wR, normalize);

opts   = optimoptions('fmincon','Display','off','Algorithm','sqp');
bounds = [0.01 1; 0.05 0.5; 0.05 0.5];   % beta, gamma, alpha
min_error = inf;
best_params = [];

for trial = 1:num_trials
    lo = max(bounds(:,1),1e-6);
    ub = bounds(:,2);
    init = exp(log(lo) + rand(size(lo)).*(log(ub)-log(lo)));
    init = init';

    [params,fval] = fmincon(objective, init, [],[],[],[], ...
        bounds(:,1), bounds(:,2), [], opts);

    if fval < min_error
        min_error = fval;
        best_params = params;
    end
end

assert(~isempty(best_params), 'fmincon returned no solution in any trial.');

tol = 1e-4;
bp = best_params(:);
if any(abs(bp - bounds(:,1)) < tol | abs(bp - bounds(:,2)) < tol)
    warning('fit_SEIR_model_OOS:binding', ...
        'Window %s-%s: parameter at bound (binding constraint) - interpret with caution.', ...
        datestr(start_date,'yyyy-mm-dd'), datestr(end_date,'yyyy-mm-dd'));
end

beta_fit  = best_params(1);
gamma_fit = best_params(2);
alpha_fit = best_params(3);
R0_est    = beta_fit / gamma_fit;

% === Training-only residuals/AIC/AICc ===
[I_train_model, R_train_model] = simulate(best_params, S0, E0, I0, R0_init, N, T_train);

resI_train = infected_train(:) - I_train_model(:);
resR_train = removed_train(:)  - R_train_model(:);

RSS_I_train = sum(resI_train.^2);
RSS_R_train = sum(resR_train.^2);

assert(RSS_I_train > 0 && RSS_R_train > 0, 'RSS = 0 in training portion -> log(-Inf).');

p = 3;
k = p + 2;
n_train = 2 * T_train;
assert(n_train - k - 1 > 0, 'Training portion too small for AICc (n_train-k-1 <= 0).');

AIC_train  = T_train*log(RSS_I_train/T_train) + T_train*log(RSS_R_train/T_train) + 2*k;
AICc_train = AIC_train + 2*k*(k+1)/(n_train-k-1);

% === Full simulation over the entire window (train+test), WITHOUT re-fitting ===
[I_full_model, R_full_model] = simulate(best_params, S0, E0, I0, R0_init, N, T);

% === Out-of-sample (test) metrics ===
infected_test = infected_fit(T_train+1:T);
removed_test  = removed_fit(T_train+1:T);
resI_test = infected_test(:) - I_full_model(T_train+1:T);
resR_test = removed_test(:)  - R_full_model(T_train+1:T);

RMSE_test_I = sqrt(mean(resI_test.^2));
RMSE_test_R = sqrt(mean(resR_test.^2));
MAE_test_I  = mean(abs(resI_test));
MAE_test_R  = mean(abs(resR_test));

% === Output structure ===
results = struct( ...
    'model','SEIR', ...
    'beta',beta_fit, 'gamma',gamma_fit, 'alpha',alpha_fit, 'R0',R0_est, ...
    'T',T, 'T_train',T_train, 'T_test',T_test, ...
    'p',p, 'k',k, 'n_train',n_train, ...
    'RSS_I_train',RSS_I_train, 'RSS_R_train',RSS_R_train, ...
    'AIC_train',AIC_train, 'AICc_train',AICc_train, ...
    'RMSE_test_I',RMSE_test_I, 'RMSE_test_R',RMSE_test_R, ...
    'MAE_test_I',MAE_test_I, 'MAE_test_R',MAE_test_R, ...
    'dates',dates_fit, ...
    'infected_real',infected_fit, 'removed_real',removed_fit, ...
    'infected_model',I_full_model, 'removed_model',R_full_model);

%% === Internal functions ===
function err = compute_error(params, S0, E0, I0, R0_0, N, T, I_data, R_data, wI, wR, use_norm)
    [I_sim, R_sim] = simulate(params, S0, E0, I0, R0_0, N, T);

    I_data = I_data(:); R_data = R_data(:);
    I_sim  = I_sim(:);  R_sim  = R_sim(:);

    if use_norm
        sI = max(eps, std(I_data));
        sR = max(eps, std(R_data));
        rI = (I_data - I_sim) / sI;
        rR = (R_data - R_sim) / sR;
        err = wI * mean(rI.^2) + wR * mean(rR.^2);
    else
        err = wI * mean((I_data - I_sim).^2) + wR * mean((R_data - R_sim).^2);
    end
end

function [I_out, R_out] = simulate(params, S0, E0, I0, R0_0, N, T)
    beta  = params(1);
    gamma = params(2);
    alpha = params(3);

    S = zeros(1,T); E = zeros(1,T); I = zeros(1,T); R = zeros(1,T);
    S(1) = S0; E(1) = E0; I(1) = I0; R(1) = R0_0;

    for t = 1:T-1
        dS = -beta * S(t) * I(t) / N;
        dE =  beta * S(t) * I(t) / N - alpha * E(t);
        dI =  alpha * E(t) - gamma * I(t);
        dR =  gamma * I(t);

        S(t+1) = S(t) + dS;
        E(t+1) = E(t) + dE;
        I(t+1) = I(t) + dI;
        R(t+1) = R(t) + dR;
    end

    I_out = I(:);
    R_out = R(:);
end

end
